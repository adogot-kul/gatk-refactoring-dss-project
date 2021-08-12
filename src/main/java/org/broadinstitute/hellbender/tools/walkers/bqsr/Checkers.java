package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Map;

import static org.broadinstitute.hellbender.tools.walkers.bqsr.AnalyzeCovariates.CSV_ARG_SHORT_NAME;
import static org.broadinstitute.hellbender.tools.walkers.bqsr.AnalyzeCovariates.PDF_ARG_SHORT_NAME;

/**
 * Utility class of AnalyzeCovariates
 */

public final class Checkers {


    private Checkers() { // class can't be instantiated
    }

    /**
     * Checks whether the last-modification-time of the inputs is consistent with their relative roles.
     *
     * This routine does not thrown an exception but may output a warning message if inconsistencies are spotted.
     *
     * @param beforeFile the before report file.
     * @param afterFile  the after report file.
     */
    static void checkInputReportFileLMT(final File beforeFile, final File afterFile, final boolean ignoreLastModificationTime) {

        if (ignoreLastModificationTime  || beforeFile == null || afterFile == null) {
            return; // nothing to do here
        } else if (beforeFile.lastModified() > afterFile.lastModified()) {
            Utils.warnUser("Last modification timestamp for 'Before' and 'After'"
                    + "recalibration reports are in the wrong order. Perhaps, have they been swapped?");
        }
    }

    /**
     * Checks that at least one output was requested.
     *
     * @throw UserException if no output was requested.
     */
    static void checkOutputRequested(final File pdfFile, final File csvFile) {
        if (pdfFile == null && csvFile == null) {
            throw new UserException("you need to request at least one output:"
                    + " the intermediate csv file (-" + CSV_ARG_SHORT_NAME + " FILE)"
                    + " or the final plot file (-" + PDF_ARG_SHORT_NAME + " FILE).");
        }
    }

    /**
     * Checks the value provided to input file arguments.
     *
     * @throw UserException if there is any problem cause by or under the end user's control
     *
     * @param name command line argument short name.
     * @param value the argument value.
     */
    static void checkInputReportFile(final String name,final File value) {
        if (value == null) {
            return;
        }

        String outString = inputReport(value);
        if (!value.exists() || !value.isFile() || !value.canRead()) {
            throw new CommandLineException.BadArgumentValue(name, outString);
        }

    }

    static String inputReport(File value) {
        String tempString = "input report '" + value;
        if (!value.exists()) {
            tempString = tempString + "' does not exists or is unreachable";
        }
        else if (!value.isFile()) {
            tempString = tempString + "' is not a regular file";
        }
        else if (!value.canRead()) {
            tempString = tempString + "' cannot be read";
        }
        return tempString;
    }

    /**
     * Checks the value provided for output arguments.
     *
     * @throw UserException if there is any problem cause by or under the end user's control
     *
     * @param name command line argument short name.
     * @param value the argument value.
     */
    static void checkOutputFileOverall(final String name, final File value) {

        final File parent = checkOutputFileValue(name, value);

        if (parent == null) {
            return;
        }

        String outString = outputFileParentDirectory(parent);
        if (!parent.exists() || !parent.isDirectory() || !parent.canWrite()) {
            throw new CommandLineException.BadArgumentValue(name, outString);
        }
    }

    static String outputFileParentDirectory(File parent) {
        String tempString = "the output file parent directory '" + parent;
        if (!parent.exists()) {
            tempString = tempString + "' does not exists or is unreachable";
        }
        else if (!parent.isDirectory()) {
            tempString = tempString + "' is not a directory";
        }
        else if (!parent.canWrite()) {
            tempString = tempString + "' cannot be written";
        }
        return tempString;
    }

    static File checkOutputFileValue(final String name, final File value) {

        if (value == null) {
            return value;
        }
        if (value.exists() && !value.isFile()) {
            throw new CommandLineException.BadArgumentValue(name, "the output file location '"
                    + value + "' exists as not a file");
        }
        File parent = value.getParentFile();

        return parent;
    }

    /**
     * Generates a flatter String array representation of recalibration argument differences.
     * @param diffs the differences to represent.
     *
     * @return never <code>null</code>, an array of the same length as the size of the input <code>diffs</code>.
     */
    static String[] reportDifferencesStringArray(final Map<String, ? extends CharSequence> diffs) {
        final String[] result = new String[diffs.size()];
        int i = 0;
        for (final Map.Entry<String, ? extends CharSequence> e : diffs.entrySet()) {
            result[i++] = capitalize(e.getKey()) + ": " + e.getValue();
        }
        return result;
    }

    /**
     * Returns the input string capitalizing the first letter.
     *
     * @param str the string to capitalize
     * @return never <code>null</code>.
     */
    private static String capitalize(final String str) {
        if (str.isEmpty()) {
            return str;
        } else {
            return Character.toUpperCase(str.charAt(0)) + str.substring(1);
        }
    }

}


