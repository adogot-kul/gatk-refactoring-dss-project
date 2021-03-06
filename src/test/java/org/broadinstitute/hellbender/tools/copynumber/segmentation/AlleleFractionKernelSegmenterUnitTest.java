package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionKernelSegmenterUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case
    
    /**
     * Generates alternate-allele-fraction-like data (similar to zero-mean multimodal test data
     * in {@link KernelSegmenterUnitTest#dataKernelSegmenter()}),
     * but introduces further segments by placing data on different chromosomes.
     */
    @DataProvider(name = "dataAlleleFractionKernelSegmenter")
    public Object[][] dataAlleleFractionKernelSegmenter() {
        final int numPoints = 10000;
        final double noiseLevel = 0.001;
        final double homFraction = 0.025;   //low hom fraction minimizes uncertainty in the changepoints coming from runs of adjacent homs near the changepoints

        final Random rng = new Random(RANDOM_SEED);
        final List<Double> minorAlleleFractions = Arrays.asList(0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25, 0.45, 0.05, 0.25);
        final List<Double> alternateAlleleFractions = IntStream.range(0, numPoints).boxed()
                .map(i -> rng.nextFloat() < homFraction
                        ? rng.nextBoolean()
                            ? 0. + noiseLevel * Math.abs(rng.nextGaussian())                                //hom ref
                            : 1. - noiseLevel * Math.abs(rng.nextGaussian())                                //hom alt
                        : rng.nextBoolean()
                            ? minorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian()          //het alt minor
                            : 1. - minorAlleleFractions.get(i / 1000) + noiseLevel * rng.nextGaussian())    //het ref minor
                .map(f -> Math.min(Math.max(f, 0.), 1.))
                .collect(Collectors.toList());             //changepoints at 999, 1999, 2999, 3999, 4999, 5999, 6999, 7999, 8999

        final List<SimpleInterval> intervals = IntStream.range(0, numPoints).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 2500 + 1),  //start a new chromosome every 2500 points, which adds additional changepoints
                        (i % 2500) + 1,
                        (i % 2500) + 1))
                .collect(Collectors.toList());

        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(intervals.stream()
                        .map(SimpleInterval::getContig)
                        .distinct()
                        .map(c -> new SAMSequenceRecord(c, 10000))
                        .collect(Collectors.toList())));

        final int globalDepth = 100;
        final List<AllelicCount> allelicCountsList = IntStream.range(0, numPoints).boxed()
                .map(i -> new AllelicCount(
                        intervals.get(i),
                        (int) ((1 - alternateAlleleFractions.get(i)) * globalDepth),
                        (int) (alternateAlleleFractions.get(i) * globalDepth)))
                .collect(Collectors.toList());
        final AllelicCountCollection allelicCounts = new AllelicCountCollection(metadata, allelicCountsList);

        final SimpleIntervalCollection segmentsExpected =
                new SimpleIntervalCollection(
                        metadata,
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1000),
                                new SimpleInterval("1", 1001, 2000),
                                new SimpleInterval("1", 2001, 2500),
                                new SimpleInterval("2", 1, 500),
                                new SimpleInterval("2", 501, 1500),
                                new SimpleInterval("2", 1501, 2500),
                                new SimpleInterval("3", 1, 1000),
                                new SimpleInterval("3", 1001, 2000),
                                new SimpleInterval("3", 2001, 2500),
                                new SimpleInterval("4", 1, 500),
                                new SimpleInterval("4", 501, 1500),
                                new SimpleInterval("4", 1501, 2500)));

        return new Object[][]{
                {allelicCounts, segmentsExpected}
        };
    }

    @Test(dataProvider = "dataAlleleFractionKernelSegmenter")
    public void testAlleleFractionKernelSegmenter(final AllelicCountCollection allelicCounts,
                                                  final SimpleIntervalCollection segmentsExpected) {
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVariance = 0.05;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 1.;
        final double numChangepointsPenaltyLogLinearFactor = 1.;

        final SimpleIntervalCollection segments = new AlleleFractionKernelSegmenter(allelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVariance, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        Assert.assertEquals(segments, segmentsExpected);
    }
}
