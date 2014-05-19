/*
 * ListPlugins
 */
package net.maizegenetics.analysis;

import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import java.io.File;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import java.awt.Frame;
import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Terry Casstevens
 */
public class ListPlugins extends AbstractPlugin {

    private PluginParameter<Boolean> myFull = new PluginParameter.Builder<Boolean>("full", false, Boolean.class).build();

    public ListPlugins(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

        Logger.getLogger("net.maizegenetics").setLevel(Level.OFF);
        Logger.getLogger("net.maizegenetics.plugindef").setLevel(Level.INFO);
    }

    @Override
    public DataSet processData(DataSet input) {

        Set<String> classes = getTasselClasses();
        List<String> temp = new ArrayList<>();
        for (String current : classes) {
            if (isPlugin(current)) {
                if (full()) {
                    temp.add(current);
                } else {
                    temp.add(Utils.getBasename(current));
                }
            }
        }

        Collections.sort(temp);

        for (String current : temp) {
            System.out.println(current);
        }

        return null;

    }

    private static Set<String> getTasselClasses() {

        String classpath = System.getProperty("java.class.path");
        String[] paths = classpath.split(File.pathSeparator);
        String tasselPath = null;
        for (String path : paths) {
            if (path.trim().length() != 0) {
                File file = new File(path);
                if (file.exists()) {
                    tasselPath = file.getAbsolutePath();
                    if (tasselPath.endsWith("sTASSEL.jar")) {
                        break;
                    }
                }
            }

        }

        Set<String> classes = new LinkedHashSet<>();
        try {
            ZipFile zFile = new ZipFile(tasselPath);
            Enumeration<? extends ZipEntry> entries = zFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory()) {
                    String name = entry.getName().replace(File.separator, ".");
                    if ((name.endsWith(".class")) && (!name.contains("$"))) {
                        name = name.substring(0, name.lastIndexOf(".class"));
                        classes.add(name);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return classes;
    }

    private static boolean isPlugin(String className) {
        try {
            Class currentMatch = Class.forName(className);
            if (!currentMatch.isInterface() && Plugin.class.isAssignableFrom(currentMatch)) {
                return true;
            } else {
                return false;
            }
        } catch (Exception ex) {
            return false;
        }
    }

    public boolean full() {
        return myFull.value();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "List Plugins";
    }

    @Override
    public String getToolTipText() {
        return "List Plugins";
    }
}
