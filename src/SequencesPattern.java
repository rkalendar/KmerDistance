import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public final class SequencesPattern {

    public SequencesPattern() {
        pt = new HashMap<>();
        pt.put("aaatt", 0);
        pt.put("aagtt", 1);
        pt.put("ataat", 2);
        pt.put("atgat", 3);
        pt.put("acagt", 4);
        pt.put("acggt", 5);
        pt.put("agact", 6);
        pt.put("aggct", 7);
        pt.put("ttaaa", 8);
        pt.put("ttgaa", 9);
        pt.put("taata", 10);
        pt.put("tagta", 11);
        pt.put("tgaca", 12);
        pt.put("tggca", 13);
        pt.put("tcaga", 14);
        pt.put("tcgga", 15);
        pt.put("ccagg", 16);
        pt.put("ccggg", 17);
        pt.put("caatg", 18);
        pt.put("cagtg", 19);
        pt.put("ctaag", 20);
        pt.put("ctgag", 21);
        pt.put("cgacg", 22);
        pt.put("cggcg", 23);
        pt.put("ggacc", 24);
        pt.put("gggcc", 25);
        pt.put("gaatc", 26);
        pt.put("gagtc", 27);
        pt.put("gtaac", 28);
        pt.put("gtgac", 29);
        pt.put("gcagc", 30);
        pt.put("gcggc", 31);
        pt.put("aaatc", 32);
        pt.put("aagtc", 33);
        pt.put("ccagt", 34);
        pt.put("acgga", 35);
        pt.put("agaca", 36);
        pt.put("atatc", 37);
        pt.put("atgct", 38);
        pt.put("ttacc", 39);
        pt.put("cagac", 40);
        pt.put("ccata", 41);
        pt.put("ccgat", 42);
        pt.put("cgaac", 43);
        pt.put("cgaat", 44);
        pt.put("ctaac", 45);
        pt.put("ctttc", 46);
        pt.put("gaccg", 47);
        pt.put("gaccc", 48);
        pt.put("gccat", 49);
        pt.put("gcgca", 50);
        pt.put("ggtta", 51);
        pt.put("ggcat", 52);
        pt.put("gtatt", 53);
        pt.put("gtgga", 54);
        pt.put("tacac", 55);
        pt.put("taggc", 56);
        pt.put("tcaaa", 57);
        pt.put("tcgcc", 58);
        pt.put("tggat", 59);
        pt.put("tgtta", 60);
        pt.put("ttggt", 61);
    }

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
    }

    public void SetFolder(String[] tfiles, String filePath) throws IOException {
        reportfile = filePath + File.separator + "result.xls";

        seq = new String[0];
        sname = new String[0];

        for (String infile : tfiles) {
            try {
                System.out.println("Target file: " + infile);
                byte[] binaryArray = Files.readAllBytes(Paths.get(infile));
                ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
                if (rf.getNseq() == 0) {
                    System.out.println("There is no sequence(s).");
                    System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                } else {
                    String[] extSeq = new String[seq.length + rf.getNseq()];
                    String[] extNam = new String[sname.length + rf.getNseq()];
                    System.arraycopy(seq, 0, extSeq, 0, seq.length);
                    System.arraycopy(sname, 0, extNam, 0, sname.length);
                    String[] s = rf.getSequences();
                    String[] n = rf.getNames();
                    System.arraycopy(s, 0, extSeq, seq.length, s.length);
                    System.arraycopy(n, 0, extNam, sname.length, n.length);
                    seq = extSeq;
                    sname = extNam;
                    if (rf.getNseq() > 1) {
                        System.out.println("Target FASTA sequences = " + rf.getNseq());
                    }
                }
            } catch (IOException e) {
                System.err.println("Failed to open file: " + infile);
            }
        }
    }

    public void RunPattern() throws IOException {
        StringBuilder sr = new StringBuilder();

        int nseq = seq.length;
        int nkmers = pt.size();
        int kmer = 5;
        int k = 0;

        String[] keys = new String[nkmers];
        for (Map.Entry<String, Integer> entry : pt.entrySet()) {
            String key = entry.getKey();
            kmer = key.length();
            Integer value = entry.getValue();
            keys[value] = key;
        }

        sr.append("ID\t");
        for (int n1 = 0; n1 < nkmers; n1++) {
            for (int n2 = n1 + 1; n2 < nkmers; n2++) {
                sr.append(keys[n1]).append("_").append(keys[n2]).append("\t");
            }
        }
        sr.append("\n");

        for (int n1 = 0; n1 < nkmers; n1++) {
            for (int n2 = n1 + 1; n2 < nkmers; n2++) {
                k++;
            }
        }

        int[][] m2 = new int[nseq][nkmers];
        int[][] m3 = new int[nseq][k];

        for (int j = 0; j < nseq; j++) {
            String r = dna.ComplementDNA2(seq[j]);
            for (int i = 0; i < seq[j].length() - kmer + 1; i++) {
                String s = seq[j].substring(i, i + kmer);
                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }
                s = r.substring(i, i + kmer);
                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }
            }

            int h = -1;
            for (int n1 = 0; n1 < nkmers; n1++) {
                int v1 = m2[j][n1];
                for (int n2 = n1 + 1; n2 < nkmers; n2++) {
                    int v2 = m2[j][n2];
                    m3[j][++h] = 0;
                    if (v1 > 0 && v2 > 0) {
                        double d;
                        if (v1 < v2) {
                            d = (1000 * v1) / v2;
                        } else {
                            d = (1000 * v2) / v1;
                        }
                        m3[j][h] = (int) d;
                    }
                }
            }
            sr.append(sname[j]).append("\t");
            for (int n2 = 0; n2 < k; n2++) {
                sr.append(m3[j][n2]).append("\t");
            }
            sr.append("\n");
        }

        int[][] mf = new int[nseq][nseq];

        for (int j = 0; j < nseq - 1; j++) {
            for (int i = j + 1; i < nseq; i++) {
                int v = 0;
                for (int n = 0; n < k; n++) {
                    if (m3[j][n] == m3[i][n]) {
                        v++;
                    } else {
                        if (m3[j][n] > m3[i][n]) {
                            if (m3[j][n] <= (dif * m3[i][n])) {
                                v++;
                            }
                        } else {
                            if ((m3[j][n] * dif) >= m3[i][n]) {
                                v++;
                            }
                        }
                    }
                }
                mf[j][i] = ((100 * v) / k);
                mf[i][j] = mf[j][i];
            }
        }

        sr.append("\n \t");
        for (int j = 0; j < nseq; j++) {
            sr.append(sname[j]).append("\t");
        }
        sr.append("\n");
        for (int i = 0; i < nseq - 1; i++) {
            sr.append(sname[i]).append("\t");
            for (int n = 0; n < nseq; n++) {
                if (i == n) {
                    sr.append(100).append("\t");
                } else {
                    sr.append(mf[i][n]).append("\t");
                }
            }
            sr.append("\n");
        }
        sr.append("\n");
        try (FileWriter fileWriter = new FileWriter(reportfile); BufferedWriter bufferedWriter = new BufferedWriter(fileWriter)) {
            bufferedWriter.write(sr.toString());
        }
    }

    private final double dif = 1.4d; //Dispersion 
    private final HashMap<String, Integer> pt;
    private String reportfile;
    private String[] seq;
    private String[] sname;

}
