package games.ptsp;

import games.math.Vector2d;

/**
 * Created by IntelliJ IDEA.
 * User: sml
 * Date: 03-Feb-2005
 * Time: 15:15:25
 * To change this template use Options | File Templates.
 */
public class Agent {
    Vector2d s, ps;
    Vector2d v;

    public Agent() {
        this(0, 0);
    }

    public Agent(int x, int y) {
        s = new Vector2d(x, y);
        ps = new Vector2d(x, y);
        v = new Vector2d(0, 0);
    }

    public void update(int d) {
        ps.set(s);
        Vector2d tmp = new Vector2d();
        Vector2d acc = Run.dir[d];

        tmp.set(v);
        tmp.mul(Run.dt);
        s.add(tmp);

        tmp.set(acc);
        tmp.mul(Run.t2);
        s.add(tmp);

        // now update velocity  v = u + at
        tmp.set(acc);
        tmp.mul(Run.a);
        v.add(tmp);
    }
}
