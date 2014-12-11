package utilities;

import java.io.*;
import java.util.*;

public class File2String {
  // convenient for the string tokenizer
  public static String get(String filename) {
    try {
      File file = new File(filename);
      byte[] b = new byte[(int) file.length()];
      DataInputStream dis =
        new DataInputStream(new FileInputStream(file));
      dis.readFully(b);
      return new String(b);
    }
    catch(Exception e) {
      System.out.println(e);
      return null;
    }
  }

  public static String read(InputStream is) {
    // dont know the size - therefore
    // use a ByteArrayOutputStream
    String s = null;
    try {
      ByteArrayOutputStream bos = new ByteArrayOutputStream();
      int b;
      while ((b = is.read()) != -1)
        bos.write(b);
      return bos.toString();
    }
    catch(Exception e) {}
    return s;
  }

  public static String[][] getArray(String file) {
    try {
      FileInputStream is = new FileInputStream(file);
      String[][] sa = readArray(is);
      is.close();
      return sa;
    }
    catch(Exception e) {
      System.out.println(e);
      return null;
    }
  }

  public static String[][] readArray(InputStream is) {
    // dont know the size - therefore
    // use a ByteArrayOutputStream
    String s = null;
    try {
      DataInputStream dis = new DataInputStream(is);
      Vector v = new Vector();

      while ((s = dis.readLine()) != null) {
        StringTokenizer st = new StringTokenizer(s);
        int n = st.countTokens();
        String[] a = new String[n];
        for (int i=0; i<n; i++)  {
          a[i] = st.nextToken();
        }
        v.addElement(a);
      }
      String[][] sa = new String[v.size()][];
      for (int i=0; i<v.size(); i++)
        sa[i] = (String[]) v.elementAt(i);
      dis.close();
      return sa;
    }
    catch(Exception e) {}
    return null;
  }

  public static boolean put(String s, String filename) {
    try {
      return put(s, new File(filename));
    }
    catch(Exception e) {
      System.out.println(e);
      return false;
    }
  }

  public static boolean put(String s, File file) {
    try {
      PrintStream ps =
        new PrintStream(new FileOutputStream(file));
      ps.print(s);
      ps.close();
      return true;
    }
    catch(Exception e) {
      System.out.println(e);
      return false;
    }
  }



}

