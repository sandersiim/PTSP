package games.ptsp;

import utilities.JEasyFrame;

import javax.swing.*;
import java.awt.*;

public class MapView extends JComponent {

    public static void main(String[] args) throws Exception {
        Map map = Map.load(args[0]);
        MapView mv = new MapView(map);
        new JEasyFrame(mv, "Map View", true);
    }

    int rad = 5;
    Map map;
    Color background = Color.black;
    Color citySeen = Color.cyan;
    Color cityUnseen = Color.red;
    Color agentColor = Color.white;
    Agent agent;
    int n = 0;

    public MapView(Map map) {
        this.map = map;
        setFocusable(true);
    }

    public void paintComponent(Graphics g) {
        if (n<10) {
            g.setColor(background);
            g.fillRect(0, 0, map.width, map.height);
            n++;
        }
        for (int i = 0; i < map.cities.length; i++) {
            int x = map.cities[i][0];
            int y = map.cities[i][1];
            if (map.visited[i]) {
                g.setColor(citySeen);
            } else {
                g.setColor(cityUnseen);
            }
            g.fillOval(x - rad, y - rad, rad * 2, rad * 2);
        }
        if (agent != null) {
            g.setColor(agentColor);
            g.drawLine((int) agent.ps.x, (int) agent.ps.y, (int) agent.s.x, (int) agent.s.y);
        }
    }

    public Dimension getPreferredSize() {
        return new Dimension(map.width, map.height);
    }

}
