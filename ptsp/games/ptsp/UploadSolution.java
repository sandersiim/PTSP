package games.ptsp;

import utilities.File2String;

import java.net.URL;
import java.net.URLConnection;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: sml
 * Date: 10-Feb-2005
 * Time: 16:20:24
 * To change this template use Options | File Templates.
 */
public class UploadSolution {
    public static String server = "http://algoval.essex.ac.uk:8080/ptsp/";

    static  String web = server + "ProcessSubmission.jsp";
    public static void main(String[] args) throws Exception {
        String mapName = "Map-10-1";
        String solution = File2String.get("route583.txt");
        upload(solution, mapName, "Test");
    }

    public static void upload(String solution, String mapName, String userName) throws Exception {

        URL u = new URL(web);

        Content content = new Content();
        content.add("solution", solution);
        content.add("mapFile", mapName);
        content.add("name", userName);
        String cont = content.get();

        URLConnection uc = u.openConnection();
        uc.setDoOutput(true);
        uc.setDoInput(true);

        System.out.println( uc );

        PrintStream out = new PrintStream(uc.getOutputStream());
        // post(System.out, cont);
        post(out, cont);
        BufferedReader in = new BufferedReader(new InputStreamReader(uc.getInputStream()));
        String line;
        while ((line = in.readLine()) != null) {
            System.out.println(line);
        }
        // out.close();

    }

    public static void post(PrintStream out, String cont) {
        // simpler than imagined!!!
        out.print(cont);
    }

    public static class Content {
        StringBuffer content;

        public void add(String name, String value) {
            if (content == null) {
                content = new StringBuffer();
                content.append(name + "=" + value);
            } else {
                content.append("&" + name + "=" + value);
            }
        }
        public String get() {
            return content.toString();
        }
    }


}
