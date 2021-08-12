package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;

@CommandLineProgramProperties(
        summary = "UNSUPPORTED.  FOR EVALUATION ONLY. Merge M2 (eval) calls with MC3 (truth)",
        oneLineSummary = "UNSUPPORTED.  FOR EVALUATION ONLY. Merge M2 calls with MC",
        programGroup = VariantEvaluationProgramGroup.class
)
@ExperimentalFeature
public class MergeMutect2CallsWithMC3 extends AbstractConcordanceWalker {

    @Argument(doc = "Merged vcf.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    protected GATKPath outputVcf;

    public static final String CENTERS_KEY = "CENTERS";
    public static final String M2_CENTER_NAME = "M2";
    public static final String MC3_REF_COUNT_KEY = "NREF";
    public static final String MC3_ALT_COUNT_KEY = "NALT";

    public static final String M2_FILTERS_KEY = "M2_FILTERS";
    public static final VCFInfoHeaderLine M2_FILTERS_HEADER_LINE =
            new VCFInfoHeaderLine(M2_FILTERS_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "M2 filters applied to variant.");

    private VariantContextWriter vcfWriter;
    private String tumorSample;

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // we assume that the truth has a single alternate allele
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));

        return sameRefAllele && containsAltAllele;
    }

    @Override
    public void onTraversalStart() {
        final Set<VCFHeaderLine> headerLines = new HashSet<>(getTruthHeader().getMetaDataInSortedOrder());

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS);

        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(M2_FILTERS_HEADER_LINE);
        tumorSample = getEvalHeader().getMetaDataLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue();
        final VCFHeader mergedHeader = new VCFHeader(headerLines, Collections.singletonList(tumorSample));
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(mergedHeader);
    }

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        final ConcordanceState concordanceState = truthVersusEval.getConcordance();

        // get AD from M2 if available, otherwise from NALT and NREF MC3 info fields
        final int[] tumorAlleleCounts = truthVersusEval.hasEval() ? truthVersusEval.getEval().getGenotype(tumorSample).getAD() :
                new int[] {truthVersusEval.getTruth().getAttributeAsInt(MC3_REF_COUNT_KEY, 0),
                        truthVersusEval.getTruth().getAttributeAsInt(MC3_ALT_COUNT_KEY, 0)};

        final Genotype tumorGenotype = new GenotypeBuilder(tumorSample, truthVersusEval.getTruthIfPresentElseEval().getAlleles())
                .AD(tumorAlleleCounts)
                .make();
        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        concordanceState.addInVCF(this, truthVersusEval, genotypes);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

    protected VariantContextBuilder makeVariantContextBuilderWithM2Center(final VariantContext mc3) {
        final List<String> centers = mc3.hasAttribute(CENTERS_KEY) ?
                mc3.getAttributeAsStringList(CENTERS_KEY, "") : new ArrayList<>();
        centers.add(M2_CENTER_NAME);
        return new VariantContextBuilder(mc3).attribute(CENTERS_KEY, centers);
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    // setter's and getter's

    public VariantContextWriter getVcfWriter() {
        return vcfWriter;
    }

    public void setVcfWriter(VariantContextWriter vcfWriter) {
        this.vcfWriter = vcfWriter;
    }
    public static String getCentersKey() {
        return CENTERS_KEY;
    }

    public static String getM2CenterName() {
        return M2_CENTER_NAME;
    }

    public static String getM2FiltersKey() {
        return M2_FILTERS_KEY;
    }

}
