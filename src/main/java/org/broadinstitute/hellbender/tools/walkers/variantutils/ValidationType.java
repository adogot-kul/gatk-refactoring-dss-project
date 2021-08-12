package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;


/**
 *  Enum ValidationType of ValidateVariants
 */

interface ValidationTypeInterface {

    public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs);
}

public enum ValidationType implements ValidationTypeInterface{

    /**
     * Makes reference to all extra-strict tests listed below.
     */
    ALL{
        public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs) {
            if(validateVariants.hasReference())
            {
                if(!rsIDs.isEmpty())
                {
                    vc.extraStrictValidation(reportedRefAllele, observedRefAllele, rsIDs);
                }
                else{
                    vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
                    vc.validateAlternateAlleles();
                    vc.validateChromosomeCounts();
                }
            }
            else{
                if (rsIDs.isEmpty())
                {
                    vc.validateAlternateAlleles();
                    vc.validateChromosomeCounts();
                }
                else{
                    vc.validateAlternateAlleles();
                    vc.validateChromosomeCounts();
                    vc.validateRSIDs(rsIDs);
                }
            }
        }
    },

    /**
     * Check whether the reported reference base in the VCF is the same as the corresponding base in the
     * actual reference.
     */
    REF{
        public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs) {
            vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
        }
    },

    /**
     * Checks whether the variant IDs exists, only relevant if the user indicates a DBSNP vcf file
     */
    IDS{
        public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs) {
            if (!rsIDs.isEmpty()) {
                vc.validateRSIDs(rsIDs);
            }
        }
    },

    /**
     * Check whether all alternative alleles participate in a genotype call of at least on sample.
     */
    ALLELES{
        public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs) {
            vc.validateAlternateAlleles();
        }
    },

    /**
     * Check that the AN and AC annotations are consistent with the number of calls, alleles and then number these
     * are called across samples.
     */
    CHR_COUNTS{
        public void applyValidationTypeOnEnum(ValidateVariants validateVariants, VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs) {
            vc.validateChromosomeCounts();
        }
    };

    /**
     * Unmodifiable set of concrete validation types.
     *
     * <p>These are all types except {@link #ALL}.</p>
     */
    public static final Set<ValidationType> CONCRETE_TYPES;

    static {
        final Set<ValidationType> cts = new LinkedHashSet<>(values().length - 1);
        for (final ValidationType v : values()) {
            if (v != ALL)
                cts.add(v);
        }
        CONCRETE_TYPES = Collections.unmodifiableSet(cts);
    }

}

