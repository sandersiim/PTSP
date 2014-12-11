package games.ptsp;

import utilities.StatisticalSummary;
import utilities.File2String;
import utilities.ElapsedTimer;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.OutputStreamWriter;
import java.io.File;
import java.util.Date;
import java.util.ArrayList;

public class SaveResults {
    public static String fileRoot = "c:/data/ptsp/maps/";
    String logFile;
    static String ext = ".txt";
    static String prefix = "sub-";
    static String comment = "// ";
    String mapName;

    public static void main(String[] args) throws Exception {
        String mapName = "Map-10-1";
        String mapFile = mapFile(mapName);
        String route = File2String.get(fileRoot + mapName + "/" + "route583.txt");
        int[] solution = Run.readForceString(route);
        Map map = Map.load(mapFile);
        Results r = Run.summary(new PrintWriter(
                new OutputStreamWriter((System.out))), map, solution);
        SaveResults saver = new SaveResults(mapName);
        saver.save(route, r);
    }

    public static ArrayList getMaps() {
        File file = new File(fileRoot);
        File[] files = file.listFiles();
        String[] a = new String[files.length];
        ArrayList al = new ArrayList();
        for (int i=0; i<files.length; i++) {
            if (files[i].isDirectory()) {
                al.add(files[i].getName());
            }
        }
        return al;
    }

    public SaveResults(String mapName) {
        this.mapName = mapName;
        logFile = logFile(mapName);
    }

    public void save(String solution, Results results) throws Exception {
        // first save the file entry
        // then update the log file
        File dir = new File(fileRoot + mapName + "/");
        dir.mkdirs();
        // dir.close();
        saveSubmission(solution, results);
        saveLog(results);
    }

    private void saveSubmission(String solution, Results results) throws Exception {
        String fileName = fileRoot + mapName + "/" + prefix + System.currentTimeMillis() + ext;
        PrintWriter out = new PrintWriter(new FileWriter(fileName));
        out.println(comment + "name: \t" + results.name);
        out.println(comment + "email: \t" + results.email);
        out.println(comment + "date: \t" + new Date());
        out.println(comment + "nVisited \t" + results.nVisited);
        out.println(comment + "nVectors \t" + results.nVectors);
        out.println(comment + "nForces \t" + results.nForces);
        out.println("Solution:" + solution);
        out.flush();
        out.close();
        System.out.println("Saved the details");
    }

    private void saveLog(Results results) throws Exception {
        FileWriter fw = new FileWriter(logFile, true);
        PrintWriter pw = new PrintWriter(fw);
        results.write(pw);
        pw.flush();
        pw.close();
    }

    public static String logFile(String mapName) {
        return fileRoot + mapName + "/" + "log.txt";
    }

    public static String mapFile(String mapName) {
        return fileRoot + mapName + "/" + "map.txt";
    }
}
