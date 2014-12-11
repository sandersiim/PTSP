package games.ptsp;

import java.io.*;
import java.util.StringTokenizer;

public class AutoController {
    int[] directions;
    int state = 0;
    public AutoController(String routeFileName) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(routeFileName));
        String line = in.readLine().trim();
        int nVectors = Integer.parseInt(line);
        this.directions = new int[nVectors];
        for (int i=0; i<nVectors; i++) {
            line = in.readLine();
            StringTokenizer st = new StringTokenizer(line);
            int dir = Integer.parseInt(st.nextToken());
            this.directions[i] = dir;
        }
        in.close();
    }

    public int getDirection() {

        /*if (key == KeyEvent.VK_DOWN) {
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
        return CENTRE;*/
        if (this.state >= this.directions.length) {
            return 10;
        }
        return directions[this.state++];
    }

    
}
