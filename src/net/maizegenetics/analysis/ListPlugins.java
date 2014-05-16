/*
 * ListPlugins
 */
package net.maizegenetics.analysis;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;

import java.io.File;

import java.util.Enumeration;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;

/**
 *
 * @author Terry Casstevens
 */
public class ListPlugins extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ListPlugins.class);

    public ListPlugins(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        Set<String> classes = getTasselClasses();
        for (String current : classes) {
            if (isPlugin(current)) {
                System.out.println(current);
            }
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
