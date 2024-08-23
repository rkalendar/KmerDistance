
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

public class KmersPattern {

    public static void main(String[] args) {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();
                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    for (String nfile : filelist) {
                        try {
                            SaveResult(nfile);
                        } catch (Exception e) {
                            System.err.println("Failed to open file: " + nfile);
                        }
                    }

                } else {
                    SaveResult(infile);
                }
            }

        } else {
            System.out.println("KmersPattern (2024) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/KmersPattern\n");

        }
    }

    private static void SaveResult(String infile) {
        try {
            long startTime = System.nanoTime();
            byte[] binaryArray = Files.readAllBytes(Paths.get(infile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("Running...");
            SequencesPattern s1 = new SequencesPattern();
            s1.SetSequences(rf.getSequences(), rf.getNames());
            s1.SetFileName(infile);
            s1.RunPattern();

            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds");
        } catch (IOException e) {
            System.out.println("Incorrect file name.");
        }
    }

}
