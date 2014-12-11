package games.ptsp;

import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

public class KeyController extends KeyAdapter {
    static int noKey = -1;
    int key = noKey;
    int CENTRE = 0;
    int UP = 1;
    int RIGHT = 2;
    int DOWN = 3;
    int LEFT = 4;

    public int getDirection() {

        if (key == KeyEvent.VK_DOWN) {
            return DOWN;
        }
        if (key == KeyEvent.VK_UP) {
            return UP;
        }
        if (key == KeyEvent.VK_RIGHT) {
            return RIGHT;
        }
        if (key == KeyEvent.VK_LEFT) {
            return LEFT;
        }
        return CENTRE;
    }

    public void keyPressed(KeyEvent e) {
        // System.out.println(e);
        key = e.getKeyCode();
    }

    public void keyReleased(KeyEvent e) {
        key = noKey;
    }
}
