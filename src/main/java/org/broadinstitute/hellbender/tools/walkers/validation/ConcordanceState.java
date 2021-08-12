package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Refactoring_enum - 13/7/21 - Arthur Dogot - source: https://stackoverflow.com/questions/1396166/refactoring-and-removing-case-statements-when-circling-over-an-enum-structure/1396266#1396266
 */

interface InterfaceConcordanceState {

    public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval); // find better name
    public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval); // find better name
    public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes); // find better name
}

public enum ConcordanceState implements InterfaceConcordanceState {

    TRUE_POSITIVE("TP") {
        public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            if(truthVersusEval.getEval().isSNP()){
                evaluateInfoFieldConcordance.setSnpCount(evaluateInfoFieldConcordance.getSnpCount()+1);
            } else if (truthVersusEval.getEval().isIndel()) {
                evaluateInfoFieldConcordance.setIndelCount(evaluateInfoFieldConcordance.getIndelCount()+1);
            }
            evaluateInfoFieldConcordance.infoDifference(truthVersusEval.getEval(), truthVersusEval.getTruth());
        }
        public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            final ConcordanceState state = truthVersusEval.getConcordance();
            Utils.validateArg(state == ConcordanceState.TRUE_POSITIVE, "This is not a true positive.");
            Concordance.tryToWrite(concordance.getTruePositivesAndFalseNegativesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getTruth(), state));
            Concordance.tryToWrite(concordance.getTruePositivesAndFalsePositivesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getEval(), state));
        }
        public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes) {
            merger.getVcfWriter().add(merger.makeVariantContextBuilderWithM2Center(truthVersusEval.getTruth()).genotypes(genotypes).make());
        }
    },
    FALSE_POSITIVE("FP") {
        public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {}
        public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            final ConcordanceState state = truthVersusEval.getConcordance();
            Utils.validateArg(state == ConcordanceState.FALSE_POSITIVE, "This is not a false positive.");
            Concordance.tryToWrite(concordance.getTruePositivesAndFalsePositivesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getEval(), state));
        }
        public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes) {
            final VariantContext m2 = truthVersusEval.getEval();
            merger.getVcfWriter().add(new VariantContextBuilder(m2.getSource(), m2.getContig(), m2.getStart(), m2.getEnd(), m2.getAlleles())
                    .attribute(MergeMutect2CallsWithMC3.getCentersKey(), MergeMutect2CallsWithMC3.getM2CenterName()).genotypes(genotypes).make());
        }
    },
    FALSE_NEGATIVE("FN") {
        public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {}
        public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            final ConcordanceState state = truthVersusEval.getConcordance();
            Utils.validateArg(state == ConcordanceState.FALSE_NEGATIVE, "This is not a false negative.");
            Concordance.tryToWrite(concordance.getTruePositivesAndFalseNegativesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getTruth(), state));
        }
        public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes) {
            merger.getVcfWriter().add(new VariantContextBuilder(truthVersusEval.getTruth()).genotypes(genotypes).make());
        }
    },
    FILTERED_TRUE_NEGATIVE("FTN") {
        public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {}
        public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            final ConcordanceState state = truthVersusEval.getConcordance();
            Utils.validateArg(state == ConcordanceState.FILTERED_TRUE_NEGATIVE, "This is not a filtered true negative.");
            Concordance.tryToWrite(concordance.getFilteredTrueNegativesAndFalseNegativesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getEval(), state));
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            concordanceFilterAnalysis(concordance,truthVersusEval);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes) {}
    },
    FILTERED_FALSE_NEGATIVE("FFN") {
        public void countInTP(EvaluateInfoFieldConcordance evaluateInfoFieldConcordance, AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {}
        public void writeInVCF(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
            final ConcordanceState state = truthVersusEval.getConcordance();
            Utils.validateArg(state == ConcordanceState.FILTERED_FALSE_NEGATIVE, "This is not a filtered false negative.");
            Concordance.tryToWrite(concordance.getTruePositivesAndFalseNegativesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getTruth(), state));
            Concordance.tryToWrite(concordance.getFilteredTrueNegativesAndFalseNegativesVcfWriter(), concordance.annotateWithConcordanceState(truthVersusEval.getEval(), state));
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            concordanceFilterAnalysis(concordance,truthVersusEval);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        public void adder(MergeMutect2CallsWithMC3 merger, AbstractConcordanceWalker.TruthVersusEval truthVersusEval, final List<Genotype> genotypes) {
            final VariantContextBuilder vcb = merger.makeVariantContextBuilderWithM2Center(truthVersusEval.getTruth())
                    .attribute(MergeMutect2CallsWithMC3.getM2FiltersKey(), truthVersusEval.getEval().getFilters().stream().collect(Collectors.toList()))
                    .genotypes(genotypes);
            merger.getVcfWriter().add(vcb.make());
        }
    };
    private final String abbreviation;

    ConcordanceState(final String abbreviation) {
        this.abbreviation = abbreviation;
    }

    public String getAbbreviation() { return abbreviation; }

    void concordanceFilterAnalysis(Concordance concordance, final AbstractConcordanceWalker.TruthVersusEval truthVersusEval) {
        if (concordance.filterAnalysis != null) {
            final Set<String> filters = truthVersusEval.getEval().getFilters();
            final boolean unique = filters.size() == 1;
            filters.stream().map(concordance.filterAnalysisRecords::get).forEach(record -> concordance.updateFilterAnalysisRecord(record, this, unique));
        }
    }
}