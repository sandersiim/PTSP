package vecmat;

import utilities.File2String;

import java.awt.*;
import java.io.*;
import java.util.*;

/* So far have only implemented the things necessary for
	the MLP implementation -- these things will need to be
	added to as time goes by */

public class VecMat {
    // this provides a number of arithmetic routines to apply
    // to vectors and matrices
    // the idea is to provide lots of methods with the same
    // name that can be appied to different combinations of arguments

    public static void main(String[] args) throws Exception {
        // test reading / writing a matrix
        double[][] m = { { 0, 1, 2}, {3, 4, 5} };
        String file = "matrix.txt";
        PrintStream out = new PrintStream( new FileOutputStream ( file ) );
        print( m, out) ; // , out );
        out.close();
        String s = File2String.get( file );
        double[][] w = readMat( s, m.length, m[0].length );
        print( w );
    }

    public static double sum(double[] x) {
        double tot = 0.0;
        for (int i=0; i<x.length; i++) {
            tot += x[i];
        }
        return tot;
    }

    public static double sum(int[] x) {
        double tot = 0.0;
        for (int i=0; i<x.length; i++) {
            tot += x[i];
        }
        return tot;
    }

    public static void invertDiagonal(double[][] m, double[][] inv) {
        // puts inverted version of m in inv
        // assumes that the matrix is diagonal
        // (and square)
        for (int i = 0; i < m.length; i++)
            inv[i][i] = 1.0 / m[i][i];
    }

    public static double diagonalDet(double[][] m) {
        // puts inverted version of m in inv
        // assumes that the matrix is diagonal
        // (and square)
        double prod = m[0][0];
        for (int i = 1; i < m.length; i++)
            prod *= m[i][i];
        return prod;
    }

    private static void vvvFast(float[] v1, float[] v2, float[] v3) {
        for (int i = 0; i < v1.length; i++)
            v3[i] = v1[i] * v2[i];
    }

    private static void vvvSafe(float[] v1, float[] v2, float[] v3) {
        int MaxLength = Math.max(v1.length, Math.max(v2.length, v3.length));
        for (int i = 0; i < MaxLength; i++)
            v3[i % v3.length] = v1[i % v1.length] * v2[i % v2.length];
    }

    public static void append(double[] v1, double[] v2, double[] v3) {
        if ((v1.length + v2.length) != v3.length) {
            System.out.println("Error in append VecMat.append");
            return;
        }
        int index = 0;
        for (int i = 0; i < v1.length; i++)
            v3[index++] = i;
        for (int i = 0; i < v2.length; i++)
            v3[index++] = v2[i];
    }

    public static double min(double[] v) {
        double m = v[0];
        for (int i = 0; i < v.length; i++)
            m = Math.min(m, v[i]);
        return m;
    }

    public static double max(double[] v) {
        double m = v[0];
        for (int i = 0; i < v.length; i++)
            m = Math.max(m, v[i]);
        return m;
    }

    public static int max(int[] v) {
        int m = v[0];
        for (int i = 0; i < v.length; i++)
            m = Math.max(m, v[i]);
        return m;
    }

    public static int argmax(float[] v) {
        int best = 0;
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; i++)
            if (v[i] > max) {
                best = i;
                max = v[i];
            }
        return best;
    }

    public static int argmax(double[] v) {
        int best = 0;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; i++)
            if (v[i] > max) {
                best = i;
                max = v[i];
            }
        return best;
    }

    public static int argmax(int[] v) {
        int best = 0;
        int max = v[0];
        for (int i = 0; i < v.length; i++)
            if (v[i] > max) {
                best = i;
                max = v[i];
            }
        return best;
    }

    public static void zero(double[] v) {
        for (int i = 0; i < v.length; i++)
            v[i] = 0.0;
    }

    public static void zero(double[][] m) {
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                m[i][j] = 0.0;
    }

    public static double[] n1(int i, int n) {
        // create a vector with only the
        double[] v = new double[n];
        v[i] = 1.0;
        return v;
    }

    public static void copy(double[] v, double[] w) {
        int len = Math.min(v.length, w.length);
        for (int i = 0; i < len; i++)
            w[i] = v[i];
    }

    public static void copy(double[][] v, double[][] w) {
        for (int i = 0; i < v.length; i++)
            for (int j = 0; j < v[0].length; j++)
                w[i][j] = v[i][j];
    }

    public static double[][] id(int n) {
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i == j)
                    m[i][j] = 1.0;
                else
                    m[i][j] = 0.0;
        return m;
    }

    public static void print(float[] v) {
        for (int i = 0; i < v.length; i++)
            System.out.print(v[i] + " : ");
        System.out.println("");
    }

    public static void print(double[] v) {
        for (int i = 0; i < v.length; i++)
            System.out.print(v[i] + " : ");
        System.out.println("");
    }

    public static void print(int[] v) {
        for (int i = 0; i < v.length; i++)
            System.out.print(v[i] + " : ");
        System.out.println("");
    }

    public static void print(float[][] v) {
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++)
                System.out.print(v[i][j] + " : ");
            System.out.println("");
        }
    }

    public static void print(double[][] v) {
        print( v , System.out );
    }

    public static void print(double[][] v, PrintStream out) {
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++)
                out.print((float) v[i][j] + " ");
            out.println("");
        }
    }

    public static void prod(float[] v1, float[] v2, float[] v3) {
        // there are different possibilities here;
        // if the sizes are all the same we assume a simple
        // one on one type of arrangement -- and we can dispense
        // with the modular arithmetic to ensure now out
        // of bounds errors.
        // otherwise we have to include it

        if ((v1.length == v2.length) && (v2.length == v3.length))
            vvvFast(v1, v2, v3);
        else
            vvvSafe(v1, v2, v3);
    }

    public static void add(float[] v1, float[] v2, float[] v3) {
        for (int i = 0; i < v1.length; i++)
            v3[i] = v1[i] + v2[i];
    }

    public static void prod(float[] v1, float[] v2, float[][] m) {
        for (int i = 0; i < v1.length; i++)
            for (int j = 0; j < v2.length; j++)
                m[i][j] = v1[i] * v2[j];
    }

    public static void prod(double[] v1, double[] v2, double[][] m) {
        for (int i = 0; i < v1.length; i++)
            for (int j = 0; j < v2.length; j++)
                m[i][j] = v1[i] * v2[j];
    }

    public static void sigmoid(float[] v1, float[] v2) {
        for (int i = 0; i < v1.length; i++)
            v2[i] = 1 / (1 + (float) Math.exp(-v1[i]));
    }

    public static void sigmoid(double[] v1, double[] v2) {
        for (int i = 0; i < v1.length; i++)
            v2[i] = 1 / (1 + Math.exp(-v1[i]));
    }

    public static double sigmoid(double x) {
        return 1.0 / (1.0 + Math.exp(-x));
    }

    public static void randomise(float[][] m) {
        // Random r = new Random();
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                m[i][j] = (float) (Math.random() * 3.0 - 1.5); // r.nextFloat();
    }

    public static void randomise(double[][] m) {
        // Random r = new Random();
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                m[i][j] = (Math.random() * 3.0 - 1.5); // r.nextFloat();
    }

    public static float applyPlus(float[] v) {
        float tot = 0;
        for (int i = 0; i < v.length; i++)
            tot += v[i];
        return tot;
    }

    public static double applyPlus(double[] v) {
        double tot = 0;
        for (int i = 0; i < v.length; i++)
            tot += v[i];
        return tot;
    }

    public static void prod(float[] v1, float[][] m, float[] v2) {
        float tot;
        if ((v1.length == m.length) && // m.length is no. of cols
                (v2.length == m[0].length)) {
            for (int i = 0; i < v2.length; i++) {
                tot = 0;
                for (int j = 0; j < v1.length; j++)
                    tot += v1[j] * m[j][i];
                v2[i] = tot;
            }
        } else
            System.out.println("Mismatch in vmv");
    }

    public static void prod(double[] v1, double[][] m, double[] v2) {
        double tot;
        if ((v1.length == m.length) && // m.length is no. of cols
                (v2.length == m[0].length)) {
            for (int i = 0; i < v2.length; i++) {
                tot = 0;
                for (int j = 0; j < v1.length; j++)
                    tot += v1[j] * m[j][i];
                v2[i] = tot;
            }
        } else
            System.out.println("Mismatch in vmv");
    }

    public static void prod(float[][] m, float[] v1, float[] v2) {
        float tot;
        if ((v1.length == m[0].length) &&
                (v2.length == m.length)) {
            for (int i = 0; i < v2.length; i++) {
                tot = 0;
                for (int j = 0; j < v1.length; j++)
                    tot += v1[j] * m[i][j];
                v2[i] = tot;
            }
        } else
            System.out.println("Mismatch in vmv");
    }

    public static double[][] readMat(String data, int cols, int rows) {
        double[][] m = new double[cols][rows];
        StringTokenizer st = new StringTokenizer( data );
        for (int i=0; i<cols; i++) {
            for(int j=0; j<rows; j++) {
                double x = new Double( st.nextToken() ).doubleValue();
                m[i][j] = x;
            }
        }
        return m;
    }

    public static void prod(double[][] m, double[] v1, double[] v2) {
        double tot;
        if ((v1.length == m[0].length) &&
                (v2.length == m.length)) {
            for (int i = 0; i < v2.length; i++) {
                tot = 0;
                for (int j = 0; j < v1.length; j++)
                    tot += v1[j] * m[i][j];
                v2[i] = tot;
            }
        } else
            System.out.println("Mismatch in vmv");
    }

    public static void add(float[][] m1, float[][] m2, float[][] m3) {
        for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1[0].length; j++)
                m3[i][j] = m1[i][j] + m2[i][j];
    }

    public static void add(double[][] m1, double[][] m2, double[][] m3) {
        for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1[0].length; j++)
                m3[i][j] = m1[i][j] + m2[i][j];
    }

    public static void prod(double[][] m1, double[][] m2, double[][] m3) {
        // no error checking - assumes all matrices are square and the same size
        int n = m1.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double tot = 0.0;
                for (int k = 0; k < n; k++)
                    tot += m1[i][k] * m2[k][j];
                m3[i][j] = tot;
            }
    }

    public static void wadd(double[][] m1, double[][] m2, double[][] m3, double w) {
        // weighted addition
        for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1[0].length; j++)
                m3[i][j] = m1[i][j] + w * m2[i][j];
    }

    public static void wadd(double[] m1, double[] m2, double[] m3, double w) {
        // weighted addition
        for (int i = 0; i < m1.length; i++)
            m3[i] = m1[i] + w * m2[i];
    }

    public static void scale(double[] v, double w) {
        for (int i = 0; i < v.length; i++)
            v[i] = v[i] * w;
    }

    public static void add(double[] m1, double[] m2, double[] m3) {
        for (int i = 0; i < m1.length; i++)
            m3[i] = m1[i] + m2[i];
    }

    public static void prod(double[] m1, double[] m2, double[] m3) {
        for (int i = 0; i < m1.length; i++)
            m3[i] = m1[i] * m2[i];
    }

    private static void safeSub(float[] v1, float[] v2, float[] v3) {
        int MaxLength = Math.max(v1.length, Math.max(v2.length, v3.length));
        for (int i = 0; i < MaxLength; i++)
            v3[i % v3.length] = v1[i % v1.length] - v2[i % v2.length];
    }

    public static void sub(float[] v1, float[] v2, float[] v3) {
        // do a safe and slow implementation for now
        // nah! assume its ok instead!!! safeSub(v1, v2, v3);
        for (int i = 0; i < v1.length; i++)
            v3[i] = v1[i] - v2[i];
    }

    public static void sub(double[] v1, double[] v2, double[] v3) {
        for (int i = 0; i < v1.length; i++)
            v3[i] = v1[i] - v2[i];
    }

    public static void sub(float s, float[] v2, float[] v3) {
        // subtract each element from a scalar s
        for (int i = 0; i < v2.length; i++)
            v3[i] = s - v2[i];
    }

    public static void sub(double s, double[] v2, double[] v3) {
        // subtract each element from a scalar s
        for (int i = 0; i < v2.length; i++)
            v3[i] = s - v2[i];
    }

    public static double dotprod(double[] v1, double[] v2) {
        double tot = 0.0;
        if (v1.length != v2.length)
            System.out.println(v1.length + " | " + v2.length);
        for (int i = 0; i < v1.length; i++)
            tot += v1[i] * v2[i];
        return tot;
    }

    public static double eucdist2(double[] v1, double[] v2) {
        // squared euclidean distance
        double tot = 0.0;
        if (v1.length != v2.length)
            System.out.println(v1.length + " | " + v2.length);
        for (int i = 0; i < v1.length; i++)
            tot += sqr(v1[i] - v2[i]);
        return tot;
    }

    public static double eucdist(double[] v1, double[] v2) {
        double tot = 0.0;
        return Math.sqrt(eucdist2(v1, v2));
    }

    /** makes the sums of the squares equal to 1.0 */
    public static void normalise(double[] v) {
        double mag = 0.0;
        for (int i = 0; i < v.length; i++)
            mag += sqr(v[i]);
        mag = Math.sqrt(mag);
        for (int i = 0; i < v.length; i++)
            v[i] = v[i] / mag;
    }

    /** makes the sums of the abs values equal to 1.0 */
    public static void normaliseAbs(double[] v) {
        double mag = 0.0;
        for (int i = 0; i < v.length; i++) {
            mag += Math.abs(v[i]);
        }
        for (int i = 0; i < v.length; i++) {
            v[i] = v[i] / mag;
        }
    }

    public static double sqr(double x) {
        return x * x;
    }
}
