package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;

/**
 * Sub class of ValidateVariants, extending the syper class in case VALIDATE_GVCF == true
 */
public class GVFCValidateVariants extends ValidateVariants {

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        super.apply(vc,readsContext, ref, featureContext);
        // if (VALIDATE_GVCF) {
            final SimpleInterval refInterval = ref.getInterval();
            // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
            // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
            // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
            // small as possible while still containing the same bases.
            final int start = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    previousInterval.getStart() : refInterval.getStart();
            final int end = (previousInterval != null && previousInterval.overlapsWithMargin(refInterval, 1)) ?
                    Math.max(previousInterval.getEnd(), vc.getEnd()) : vc.getEnd();
            final GenomeLoc possiblyMergedGenomeLoc = genomeLocSortedSet.getGenomeLocParser().createGenomeLoc(refInterval.getContig(), start, end);
            genomeLocSortedSet.add(possiblyMergedGenomeLoc, true);

            previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
            previousStart = vc.getStart();
            validateGVCFVariant(vc);
        //}
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        // if (VALIDATE_GVCF) {
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

            if (seqDictionary == null)
                throw new UserException("Validating a GVCF requires a sequence dictionary but no dictionary was able to be constructed from your input.");

            genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        // }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        // if (VALIDATE_GVCF) {
            final GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

            if (intervalArgumentCollection.intervalsSpecified()){
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
            } else {
                intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(seqDictionary);
            }

            final GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(genomeLocSortedSet);
            if (uncoveredIntervals.coveredSize() > 0) {
                final UserException e = new UserException("A GVCF must cover the entire region. Found " + uncoveredIntervals.coveredSize() +
                        " loci with no VariantContext covering it. The first uncovered segment is:" +
                        uncoveredIntervals.iterator().next());
                throwOrWarn(e);
            }
       //  }
        return null;
    }
}
