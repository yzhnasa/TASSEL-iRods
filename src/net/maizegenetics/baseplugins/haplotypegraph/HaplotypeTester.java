package net.maizegenetics.baseplugins.haplotypegraph;

import javax.swing.*;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 6, 2004
 * Time: 8:44:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class HaplotypeTester
{
    public static void main(String[] args) throws Exception
    {
        UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        HaplotypeFrame hf = new HaplotypeFrame( new File("C:\\IdeaProjects\\Tassel\\haplotyeTest.txt") );
        hf.run();
    }
}
