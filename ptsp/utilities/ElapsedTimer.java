package utilities;

public class ElapsedTimer {
    // allows for easy reporting of elapsed time
    long oldTime;

    public ElapsedTimer() {
        oldTime = System.currentTimeMillis();
    }

    public long elapsed() {
        return System.currentTimeMillis() - oldTime;
    }

    public void reset() {
        oldTime = System.currentTimeMillis();
    }

    public String toString() {
        return elapsed() + " ms elapsed";
    }

    public static void main(String[] args) {

        ElapsedTimer t = new ElapsedTimer();

        System.out.println("ms elasped: " + t.elapsed());
    }

}


