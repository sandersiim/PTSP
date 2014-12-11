package games.ptsp;

import games.math.Vector2d;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Random;
import java.net.URL;
import java.net.URLConnection;

public class Map {
    static int width = 320;
    static int height = 240;
    // static int nCities = 20;
    static int seed = 1;
    static int empty = 0;
    static int city = 1;
    static Random r = new Random();
    boolean[] visited;
    int[][] cities;
    public int nCities;

    public static void main(String[] args) throws Exception {
        // make a map and save it to a file
        int nCities = 10;
        Map m = new Map(nCities);
        m.save("cities.txt");
    }

    // 2d map of all cities
    // list of x-y city coords

    public Map(int nCities) {
        this.nCities = nCities;
        cities = new int[nCities][];
        visited = new boolean[nCities];
        randomise();
    }

    public void save(String filename) throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(filename));
        out.close();
    }

    public void save(PrintWriter out) {
        out.println(nCities);
        for (int i=0; i<nCities; i++) {
            out.println(cities[i][0] + " " + cities[i][1]);
        }
    }

    public static Map load(String filename) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(filename));
        String line = in.readLine().trim();
        int nCities = Integer.parseInt(line);
        Map map = new Map(nCities);
        for (int i=0; i<nCities; i++) {
            line = in.readLine();
            StringTokenizer st = new StringTokenizer(line);
            int x = Integer.parseInt(st.nextToken());
            int y = Integer.parseInt(st.nextToken());
            map.cities[i] = new int[]{x, y};
        }
        in.close();
        return map;
    }

    public static Map load(BufferedReader in) throws Exception {
        String line = in.readLine().trim();
        int nCities = Integer.parseInt(line);
        Map map = new Map(nCities);
        for (int i=0; i<nCities; i++) {
            line = in.readLine();
            StringTokenizer st = new StringTokenizer(line);
            int x = Integer.parseInt(st.nextToken());
            int y = Integer.parseInt(st.nextToken());
            map.cities[i] = new int[]{x, y};
        }
        return map;
    }

    public static Map load(URL u) throws Exception {
        URLConnection uc = u.openConnection();
        InputStreamReader isr = new InputStreamReader(uc.getInputStream());
        BufferedReader in = new BufferedReader(isr);
        // get rid of blank line at beginning...
        in.readLine();
        Map map = load(in);
        in.close();
        return map;
    }

    public boolean tryVisit(double x, double y) {
        for (int i=0; i<nCities; i++) {
            if (!visited[i] && cities[i] != null) {
                if (vicinity(cities[i], x, y)) {
                    // checks then makes the visit
                    visited[i] = true;
                    return true;
                }
            }
        }
        return false;
    }

    public boolean checkVisit(double x, double y) {
        for (int i=0; i<nCities; i++) {
            if (!visited[i] && cities[i] != null) {
                if (vicinity(cities[i], x, y)) {
                    // visited[i] = true;
                    return true;
                }
            }
        }
        return false;
    }

    public boolean vicinity(int[] city, double x, double y) {
        double xd = city[0] - x;
        double yd = city[1] - y;
        return xd * xd + yd * yd < Run.rad * Run.rad;
    }

    public void reset() {
        visited = new boolean[nCities];
    }

    public void randomise() {
        int nMade = 0;
        while (nMade < nCities) {
            int x = r.nextInt(width);
            int y = r.nextInt(height);
            if (!checkVisit(x, y)) {
                // place a new city here
                cities[nMade] = new int[]{x, y};
                nMade++;
            }
        }
    }
}
