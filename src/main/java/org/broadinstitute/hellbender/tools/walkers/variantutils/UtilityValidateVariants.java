package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.*;

/**
 * Utility class of ValidateVariants
 */
public final class UtilityValidateVariants {

    private UtilityValidateVariants() { // class can't be instantiated
    }

    static final Logger logger = LogManager.getLogger(UtilityValidateVariants.class);

    /**
     * Given the validation type and exclusion type, calculate the final set of type to validate.
     * @param excludeTypes types to exclude.
     *
     * @return the final set of type to validate. May be empty.
     */

    /**
     * Contains final set of validation to apply.
     */


    protected static Collection<ValidationType> calculateValidationTypesToApply(ValidateVariants validateVariants) {

        final Set<ValidationType> excludeTypeSet = usesExcludeTypesTemp(validateVariants);
        if (excludeTypeSet.contains(ValidationType.ALL)) {
            if (excludeTypeSet.size() > 1) {
                logger.warn("found ALL in the --validation-type-to-exclude list together with other concrete type exclusions that are redundant");
            }
            return Collections.emptyList();
        } else {
            final Set<ValidationType> result = new LinkedHashSet<>(ValidationType.CONCRETE_TYPES);
            result.removeAll(excludeTypeSet);
            if (result.contains(ValidationType.REF) && !validateVariants.hasReference()) {
                throw new UserException.MissingReference("Validation type " + ValidationType.REF.name() + " was selected but no reference was provided.", true);
            }
            return result;
        }

    }

    private static Set<ValidationType> usesExcludeTypesTemp(ValidateVariants validateVariants) {
        //creates local, temp list so that original list provided by user doesn't get modified
        List<ValidationType> excludeTypesTemp = new ArrayList<>(validateVariants.excludeTypes);
        if (validateVariants.VALIDATE_GVCF && !excludeTypesTemp.contains(ValidationType.ALLELES)) {
            // Note: in a future version allele validation might be OK for GVCFs, if that happens
            // this will be more complicated.
            logger.warn("GVCF format is currently incompatible with allele validation. Not validating Alleles.");
            excludeTypesTemp.add(ValidationType.ALLELES);
        }
        if (excludeTypesTemp.isEmpty()) {
            return Collections.singleton(ValidationType.ALL);
        }

        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeTypesTemp);
        if (excludeTypesTemp.size() != excludeTypeSet.size()) {
            logger.warn("found repeat redundant validation types listed using the --validation-type-to-exclude argument");
        }
        return excludeTypeSet;
    }


}
