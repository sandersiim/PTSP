package games.ptsp;

import games.math.Vector2d;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: sml
 * Date: 03-Feb-2005
 * Time: 13:12:55
 * To change this template use Options | File Templates.
 */
public class Run {
    static int maxRouteLength = 50000;

    public static void main(String[] args) throws Exception {

        String mapfile = args[0];
        String forcefile = args[1];
        Map map = Map.load(mapfile);
        int[] solution = readForces(forcefile);
        if (args.length == 3) {
            verbose = true;
        }
        int nVis = nVisited(map, solution);

        System.out.println("Visited: " + nVis + " out of " + map.nCities);
        System.out.println("Valid:   " + (nVis == map.nCities));
        System.out.println("Length:  " + solution.length);
        System.out.println("nForces: " + nForces(solution));
    }

    public static Results summary(PrintWriter out, Map map, int[] solution) {
        int nVis = nVisited(map, solution);
        out.println("Visited: " + nVis + " out of " + map.nCities);
        out.println("Valid:   " + (nVis == map.nCities));
        out.println("Length:  " + solution.length);
        out.println("nForces: " + nForces(solution));
        return new Results(nVis, solution.length, nForces(solution));
    }

    public static Results getResults(Map map, int[] solution) {
        int nVis = nVisited(map, solution);
        return new Results(nVis, solution.length, nForces(solution));
    }

    public static int nForces(int[] solution) {
        int tot = 0;
        for (int i=0; i<solution.length; i++) {
            if (solution[i] != 0) {
                tot++;
            }
        }
        return tot;
    }

    public static int nVisited(Map map, int[] solution) {
        // update the position: s = ut + 0.5 at^2
        Vector2d pos = new Vector2d(Map.width/2, Map.height/2);
        Vector2d vel = new Vector2d(0, 0);
        Vector2d tmp = new Vector2d();
        int nVisited = 0;
        for (int i=0; i<solution.length; i++) {
            Vector2d acc = dir[solution[i]];
            tmp.set(vel);
            tmp.mul(dt);
            pos.add(tmp);

            tmp.set(acc);
            tmp.mul(t2);
            pos.add(tmp);

            // now update velocity  v = u + at
            tmp.set(acc);
            tmp.mul(a);
            vel.add(tmp);

            // now check to see if we've visitied another city
            if (map.tryVisit(pos.x, pos.y)) {
                nVisited++;
            };
        }
        return nVisited;
    }

    public static int[] readForceString(String forces) {
        StringTokenizer st = new StringTokenizer(forces);
        int n = Integer.parseInt(st.nextToken());
        if (n > maxRouteLength) {
            throw new RuntimeException("Max route length of " + maxRouteLength + " exceeded: " + n);
        }
        int[] a = new int[n];

        for (int i=0; i<a.length; i++) {
            a[i] = Integer.parseInt(st.nextToken());
        }

        return a;
    }



    public static int[] readForces(String filename) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(filename));
        String line = in.readLine().trim();
        int n = Integer.parseInt(line);
        int[] f = new int[n];
        for (int i=0; i<n; i++) {
            line = in.readLine().trim();
            f[i] = Integer.parseInt(line);
        }
        in.close();
        return f;
    }

    static double dt = 0.1;
    static double a = 1.0;
    static double t2 = 0.5 * dt * dt;
    static double rad = 5;
    static boolean verbose = false;

    static Vector2d[] dir = new Vector2d[]{
        new Vector2d(0, 0),
        new Vector2d(0, -a),
        new Vector2d(a, 0),
        new Vector2d(0, a),
        new Vector2d(-a, 0),
    };
}
