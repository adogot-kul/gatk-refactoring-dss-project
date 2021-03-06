package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 * Utility class that allows easy creation of destinations for the HaplotypeBAMWriters
 *
 */
public abstract class HaplotypeBAMDestination {
    private final SAMFileHeader bamOutputHeader;
    private final Optional<String> haplotypeReadGroupID;
    private final static String haplotypeSampleTag = "HC";

    /**
     * Create a new HaplotypeBAMDestination
     *
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination.
     * @param haplotypeReadGroupID read group ID used when writing haplotypes as reads. Empty if
     *                             {@link HaplotypeBAMWriter.WriterType} == NO_HAPLOTYPES.
     */
    protected HaplotypeBAMDestination(SAMFileHeader sourceHeader, final Optional<String> haplotypeReadGroupID) {
        Utils.nonNull(sourceHeader, "sourceHeader cannot be null");
        Utils.nonNull(haplotypeReadGroupID, "haplotypeReadGroupID cannot be null");
        this.haplotypeReadGroupID = haplotypeReadGroupID;

        bamOutputHeader = new SAMFileHeader();
        bamOutputHeader.setSequenceDictionary(sourceHeader.getSequenceDictionary());
        bamOutputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        readGroups.addAll(sourceHeader.getReadGroups()); // include the original read groups...

        // plus an artificial read group for the haplotypes
        if (haplotypeReadGroupID.isPresent()){
            final SAMReadGroupRecord rgRec = new SAMReadGroupRecord(haplotypeReadGroupID.get());
            rgRec.setSample(haplotypeSampleTag);
            rgRec.setSequencingCenter("BI");
            readGroups.add(rgRec);
        }

        bamOutputHeader.setReadGroups(readGroups);
        final List<SAMProgramRecord> programRecords = new ArrayList<>(sourceHeader.getProgramRecords());
        programRecords.add(new SAMProgramRecord("HaplotypeBAMWriter"));
        bamOutputHeader.setProgramRecords(programRecords);
    }

    /**
     * Write a read to the output specified by this destination.
     *
     * @param read the read to write out
     */
    public abstract void add(final GATKRead read);

    /**
     * Get the read group ID that is used by this writer when writing halpotypes as reads.
     * If {@link HaplotypeBAMWriter.WriterType} is NO_HAPLOTYPES, then returns Optional.empty()
     *
     * @return read group ID
     */
    public Optional<String> getHaplotypeReadGroupID() {
        return haplotypeReadGroupID;
    }
    /**
     * Get the sample tag that is used by this writer when writing halpotypes as reads.
     *
     * @return sample tag
     */
    public String getHaplotypeSampleTag() { return haplotypeSampleTag; }

    /**
     * Close the destination
     *
     */
    abstract void close();

    /**
     * Get the SAMFileHeader that is used for writing the output for this destination.
     * @return output SAMFileHeader
     */
    public SAMFileHeader getBAMOutputHeader() {
        return bamOutputHeader;
    }

}