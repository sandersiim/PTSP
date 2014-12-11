package games.ptsp;

import utilities.ElapsedTimer;

import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;

import games.ptsp.SaveResults;

public class Results implements Comparable {

    public static void main(String[] args) throws Exception {
        Results r1 = new Results(10, 300, 200);
        Results r2 = new Results(10, 250, 200);
        Results r3 = new Results(10, 251, 100);
        Results r4 = new Results(10, 251, 10);
        Results r5 = new Results(9, 10, 10);
        Results[] ra = new Results[]{r1, r2, r3, r4, r5};
        Arrays.sort(ra);
        for (int i=0; i<ra.length; i++) {
            System.out.println(ra[i]);
        }
    }

    // static String logFile = SaveResults.logFile;

    public int nVisited;
    public int nVectors;
    public int nForces;
    public String email;
    public String name;
    public String date;

    public static Collection getResults(String mapName) throws Exception {
        return readResults(SaveResults.logFile(mapName));
    }

    public static SortedSet readResults(String logFile) throws Exception {
        SortedSet set = new TreeSet();

        BufferedReader br = new BufferedReader(new FileReader(logFile));
        String line;
        while ((line = br.readLine()) != null) {
            try {
                Results r = new Results(line);
                set.add(r);
            } catch (Exception e) {
            }
        }
        return set;
    }

    public Results(int nVisited, int nVectors, int nForces) {
        this.nVisited = nVisited;
        this.nVectors = nVectors;
        this.nForces = nForces;
    }

    public Results(int nVisited, int nVectors, int nForces, String email, String name) {
        this.nVisited = nVisited;
        this.nVectors = nVectors;
        this.nForces = nForces;
        this.email = email;
        this.name = name;
    }

    public Results(String line) {
        StringTokenizer st = new StringTokenizer(line, "\t");
        nVisited = Integer.parseInt(st.nextToken().trim());
        nVectors = Integer.parseInt(st.nextToken().trim());
        nForces = Integer.parseInt(st.nextToken().trim());
        email = st.nextToken().trim();
        name = st.nextToken().trim();
        date = st.nextToken().trim();
    }

    public void write(PrintWriter out) {
        out.println(nVisited + "\t " + nVectors + "\t " + nForces + "\t"
                + email + "\t " + name + "\t " + new Date() );
    }

    public int compareTo(Object o) {
        Results r = (Results) o;
        int c1 = -sgn(nVisited - r.nVisited);
        if (c1 != 0) {
            return c1;
        }
        int c2 = sgn(nVectors - r.nVectors);
        if (c2 != 0) {
            return c2;
        }
        return sgn(nForces - r.nForces);
    }

    public static int sgn(double x) {
        if (x == 0) return 0;
        if (x > 0) {
            return 1;
        } else {
            return -1;
        }
    }

    public String toString() {
        return nVisited + "\t " + nVectors + "\t " + nForces;
    }
}
