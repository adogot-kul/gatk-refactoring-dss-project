package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.variantutils.UtilityValidateVariants.calculateValidationTypesToApply;


/**
 * Validate a VCF file with a strict set of criteria
 *
 * <p> This tool is designed to validate the adherence of a file to VCF format. The tool will validate .g.vcf GVCF
 * format files as well. For VCF specifications, see
 * <a href='https://samtools.github.io/hts-specs/'>https://samtools.github.io/hts-specs/</a>.
 * Besides standard adherence to the VCF specification, this tool performs additional strict validations to ensure
 * that the information contained within the file is correctly encoded. These include:
 * </p>
 *
 * <ul>
 *   <li><b>REF</b> - correctness of the reference base(s)</li>
 *   <li><b>CHR_COUNTS</b> - accuracy of AC and AN values</li>
 *   <li><b>IDS</b> - tests against rsIDs when a dbSNP file is provided (requires a dbsnp VCF provided via `--dbsnp`).</li>
 *   <li><b>ALLELES</b> - that all alternate alleles are present in at least one sample</li>
 * </ul>
 *
 * <p>
 *     By default the tool applies all the strict validations unless you indicate which one should be
 *     excluded using `--validation-type-to-exclude`. You can exclude as many types as you want. Furthermore, you
 *     can exclude all strict validations with the special code `ALL`. In this case the tool will only test for
 *     adherence to the VCF specification.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A VCF file to validate.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Minimally validate a file for adherence to VCF format:</h4>
 * gatk ValidateVariants \
 *     -V cohort.vcf.gz
 *
 * <h4>Validate a GVCF for adherence to VCF format, including REF allele match:</h4>
 * gatk ValidateVariants \
 *     -V sample.g.vcf.gz \
 *     -R reference.fasta
 *     -gvcf
 *
 * <h4>To perform VCF format and all strict validations: </h4>
 * <pre>
 * gatk ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>To perform only VCF format tests:</h4>
 * <pre>
 * gatk ValidateVariants
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validation-type-to-exclude ALL
 * </pre>
 *
 * <h4>To perform all validations except the strict `ALLELE` validation:</h4>
 * <pre>
 * gatk ValidateVariants \
 *   -R ref.fasta \
 *   -V input.vcf \
 *   --validation-type-to-exclude ALLELES \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Validates a VCF file with an extra strict set of criteria.",
        oneLineSummary = "Validate VCF",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class ValidateVariants extends VariantWalker {


    public static final String GVCF_VALIDATE = "validate-GVCF";
    public static final String DO_NOT_VALIDATE_FILTERED_RECORDS = "do-not-validate-filtered-records";


    @ArgumentCollection
    DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Argument(fullName = "validation-type-to-exclude",
            shortName = "Xtype",
            doc = "which validation type to exclude from a full strict validation",
            optional = true)
    List<ValidationType> excludeTypes = new ArrayList<>();

    /**
     * By default, even filtered records are validated.
     */
    @Argument(fullName = DO_NOT_VALIDATE_FILTERED_RECORDS,
            shortName = "do-not-validate-filtered-records",
            doc = "skip validation on filtered records",
            optional = true,
            mutex = GVCF_VALIDATE)
    Boolean DO_NOT_VALIDATE_FILTERED = false;

    @Argument(fullName = "warn-on-errors",
            shortName = "warn-on-errors",
            doc = "just emit warnings on errors instead of terminating the run at the first instance",
            optional = true)
    Boolean WARN_ON_ERROR = false;

    /**
     *  Validate this file as a gvcf. In particular, every variant must include a <NON_REF> allele, and that
     *  every base in the territory under consideration is covered by a variant (or a reference block).
     *  If you specifed intervals (using -L or -XL) to restrict analysis to a subset of genomic regions,
     *  those intervals will need to be covered in a valid gvcf.
     */
    @Argument(fullName = GVCF_VALIDATE,
            shortName = "gvcf",
            doc = "Validate this file as a GVCF",
            optional = true,
            mutex = DO_NOT_VALIDATE_FILTERED_RECORDS)
    Boolean VALIDATE_GVCF = false;


    protected GenomeLocSortedSet genomeLocSortedSet;

    // information to keep track of when validating a GVCF
    protected SimpleInterval previousInterval;
    protected int previousStart = -1;
    private String previousContig = null;

    protected Collection<ValidationType> validationTypes;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if (DO_NOT_VALIDATE_FILTERED && vc.isFiltered()) {
            return;
        }
        validateVariantsOrder(vc);
        loopThroughCollectionForValidation(vc, ref, featureContext);
    }

    public void loopThroughCollectionForValidation(final VariantContext vc, final ReferenceContext ref, final FeatureContext featureContext) {

        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();

        final Allele observedRefAllele = hasReference() ? Allele.create(Arrays.copyOf(ref.getBases(), refLength)) : null;

        final Set<String> rsIDs = getRSIDs(featureContext);
        for (final ValidationType t : validationTypes) {
            try{
                applyValidationType(vc, reportedRefAllele, observedRefAllele, rsIDs, t);
            } catch (TribbleException e) {
                throwOrWarn(new UserException.FailsStrictValidation(drivingVariantFile.getRawInputString(), t, e.getMessage()));
            }
        }
    }

    /**
     *  Returns the list of RSIDs overlapping the current variant that we're walking over.
     *  If there's no RSID or if there was not dbsnp file passed in as an argument,
     *  an empty set is returned (and then no validation is performed, see applyValidationType.
     */
    private Set<String> getRSIDs(FeatureContext featureContext) {
        Set<String> rsIDs = new LinkedHashSet<>();
        for (VariantContext rsID : featureContext.getValues(dbsnp.dbsnp)) {
            rsIDs.addAll(Arrays.asList(rsID.getID().split(VCFConstants.ID_FIELD_SEPARATOR)));
        }
        return rsIDs;
    }

    private void validateVariantsOrder(final VariantContext vc) {
        // Check if the current VC belongs to the same contig as the previous one.
        // If not, reset the start position to -1.
        if (previousContig == null || !previousContig.equals(vc.getContig())) {
            previousContig = vc.getContig();
            previousStart = -1;
        }

        //if next VC refers to a previous genomic position, throw an error
        //Note that HaplotypeCaller can emit variants that start inside of a deletion on another haplotype,
        // making v2's start less than the deletion's end
        if (previousStart > -1 && vc.getStart() < previousStart) {
            final UserException e = new UserException(String.format("In a GVCF all records must ordered. Record: %s covers a position previously traversed.",
                    vc.toStringWithoutGenotypes()));
            throwOrWarn(e);
        }
    }

    protected void validateGVCFVariant(final VariantContext vc) {
        if (!vc.hasAllele(Allele.NON_REF_ALLELE)) {
            final UserException e = new UserException(String.format("In a GVCF all records must contain a %s allele. Offending record: %s",
                    Allele.NON_REF_STRING, vc.toStringWithoutGenotypes()));
            throwOrWarn(e);
        }
    }

    private void applyValidationType(VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs, ValidationType t) {
        // Note: VariantContext.validateRSIDs blows up on an empty list (but works fine with null).
        // The workaround is to not pass an empty list.

        t.applyValidationTypeOnEnum(this, vc, reportedRefAllele, observedRefAllele, rsIDs);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalStart() {
        validationTypes = calculateValidationTypesToApply( this);
        lackofArgumentsWarnings();
    }


    public void lackofArgumentsWarnings() {
        //warn user if certain requested validations cannot be done due to lack of arguments
        if(dbsnp.dbsnp == null && (validationTypes.contains(ValidationType.ALL) || validationTypes.contains(ValidationType.IDS)))
        {
            logger.warn("IDS validation cannot be done because no DBSNP file was provided");
            logger.warn("Other possible validations will still be performed");
        }
        if(!hasReference() && (validationTypes.contains(ValidationType.ALL) || validationTypes.contains(ValidationType.REF)))
        {
            logger.warn("REF validation cannot be done because no reference file was provided");
            logger.warn("Other possible validations will still be performed");
        }
    }


    public Object onTraversalSuccess() {
        return null;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    protected void throwOrWarn(UserException e) {
        if (WARN_ON_ERROR) {
            logger.warn("***** " + e.getMessage() + " *****");
        } else {
            throw e;
        }
    }

}
