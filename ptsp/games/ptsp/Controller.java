package games.ptsp;

import utilities.JEasyFrame;

import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.net.URL;

public class Controller extends Thread {
    // this one puts it all together
    // loads in a map and so on...
    Map map;
    Agent agent;
    MapView view;
    JEasyFrame frame;
    KeyController kc;
    AutoController ac;
    static int delay = 50;
    ArrayList route;
    static String saveRouteFile = "new_route.txt";
    static String server = UploadSolution.server;
    //String userName;
    String mapName;
    String routeFile;
    boolean auto = false;


    public static void main(String[] args) throws Exception {
        String mapName = args[0];
        if (args.length > 1) {
            String routeFileName = args[1];
            Controller c = new Controller(mapName, routeFileName);
        } else {
            Controller c = new Controller(mapName);
        }
    }

    public Controller(String mapName) throws Exception {
        this.mapName = mapName;
        //this.userName = userName;
        //URL u = new URL(server + "Map.jsp?map=" + mapName);
        map = Map.load(mapName);
        view = new MapView(map);
        frame = new JEasyFrame(view, "Keyboard Interface", true);
        agent = new Agent(map.width / 2, map.height / 2);
        view.agent = agent;
        kc = new KeyController();
        view.addKeyListener(kc);
        view.requestFocus();
        route = new ArrayList();
        start();
    }

    public Controller(String mapName, String routeFileName) throws Exception {
        this.mapName = mapName;
        //URL u = new URL(server + "Map.jsp?map=" + mapName);
        map = Map.load(mapName);
        view = new MapView(map);
        frame = new JEasyFrame(view, "Keyboard Interface", true);
        agent = new Agent(map.width / 2, map.height / 2);
        view.agent = agent;
        ac = new AutoController(routeFileName);
        view.requestFocus();
        route = new ArrayList();
        this.auto = true;
        start();
    }

    public void run() {
        int n = 0;
        boolean started = false;
        frame.setTitle("PTSP!");
        while (n < map.nCities) {
            try {
                int dir;
                if (this.auto) {
                    dir = ac.getDirection();
                } else {
                    dir = kc.getDirection();
                }
                if (dir != 0) {
                    started = true;
                } 
                if (dir == 10) {
                    break;
                }
                if (started) {
                    route.add(new Integer(dir));
                    agent.update(dir);
                    if (map.tryVisit(agent.s.x, agent.s.y)) {
                        n++;
                    }
                    ;
                    view.paint(view.getGraphics());
                    frame.setTitle("Visited: " + n + " in " + route.size() + " steps");
                }
                sleep(delay);
            } catch (Exception e) {
                System.out.println(e);
            }
        }
        try {
            //saveRoute(route);
            //UploadSolution.upload(routeString(route), mapName, userName);
            System.out.println("Uploaded solution");
        } catch (Exception e) {
            System.out.println(e);
        }
    }

    public String routeString(ArrayList route) {
        StringBuffer sb = new StringBuffer();
        sb.append(route.size() + "\n");
        for (int i = 0; i < route.size(); i++) {
            sb.append(route.get(i) + "\n");
        }
        return sb.toString();
    }

    public void saveRoute(ArrayList route) {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(saveRouteFile));
            out.println(route.size());
            for (int i = 0; i < route.size(); i++) {
                out.println(route.get(i));
            }
            out.close();
        } catch (Exception e) {
            System.out.println("Exception saving route: " + route);
        }
    }
}
