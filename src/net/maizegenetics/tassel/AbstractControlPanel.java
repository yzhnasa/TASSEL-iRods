package net.maizegenetics.tassel;

import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.progress.ProgressPanel;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JMenu;
import javax.swing.JPanel;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.HashMap;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
public abstract class AbstractControlPanel extends JPanel implements ActionListener {

    protected TASSELMainFrame theTASSELMainFrame;
    protected DataTreePanel theDataTreePanel;
    protected HashMap<JButton, Plugin> buttonHash = new HashMap<JButton, Plugin>();

    public AbstractControlPanel(TASSELMainFrame theQAF, DataTreePanel theDTP) {
        theTASSELMainFrame = theQAF;
        theDataTreePanel = theDTP;
        FlowLayout layout = new FlowLayout(FlowLayout.CENTER, 0, 0) {

            @Override
            public Dimension preferredLayoutSize(Container target) {

                int maxWidth = (target.getSize().width == 0) ? Integer.MAX_VALUE : target.getSize().width;

                int resultWidth = 0;
                int resultHeight = 0;

                int currentWidth = 0;
                int currentHeight = 0;
                Component[] components = target.getComponents();
                for (int i = 0; i < components.length; i++) {

                    if (components[i].isVisible()) {
                        Dimension currentSize = components[i].getPreferredSize();

                        if (currentWidth + currentSize.width > maxWidth) {
                            resultWidth = Math.max(resultWidth, currentWidth);
                            resultHeight += currentHeight;
                            currentWidth = 0;
                            currentHeight = 0;
                        }

                        currentWidth += currentSize.width;
                        currentHeight = Math.max(currentHeight, currentSize.height);
                    }

                }

                resultWidth = Math.max(resultWidth, currentWidth);
                resultHeight += currentHeight;

                return new Dimension(resultWidth, resultHeight);
            }
        };
        setLayout(layout);
    }

    public void setDataTreePanel(DataTreePanel theDataTreePanel) {
        this.theDataTreePanel = theDataTreePanel;
    }

    protected void addPlugin(Plugin theTP) {
        JButton theButton = createButton(theTP);
        theButton.addActionListener(this);
        theTP.addListener(theDataTreePanel);
        this.add(theButton);
        buttonHash.put(theButton, theTP);
    }

    /**
     * This is to add plugins which will have a JPanel displayed
     * in TASSEL's main workspace.
     */
    protected void addMainDisplayPlugin(final Plugin theTP) {

        JButton theButton = createButton(theTP);

        theButton.addActionListener(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                theTASSELMainFrame.updateMainDisplayPanel(theTP.getPanel());
            }
        });

        this.add(theButton);
        buttonHash.put(theButton, theTP);

        JMenu menu = theTP.getMenu();
        if (menu != null) {
            theTASSELMainFrame.addMenu(menu);
        }

    }

    private JButton createButton(Plugin theTP) {
        JButton theButton = new JButton(theTP.getButtonName(), theTP.getIcon());
        theButton.setFont(new Font("Dialog", 1, 12));
        theButton.setBackground(Color.white);
        theButton.setMargin(new Insets(2, 2, 2, 2));
        theButton.setToolTipText(theTP.getToolTipText());
        return theButton;
    }

    public void actionPerformed(ActionEvent e) {
        JButton theButton = (JButton) e.getSource();
        Plugin theTP = this.buttonHash.get(theButton);
        PluginEvent event = new PluginEvent(theDataTreePanel.getSelectedTasselDataSet());
        ProgressPanel progressPanel = theTASSELMainFrame.getProgressPanel();
        progressPanel.addPlugin(theTP);
        ThreadedPluginListener thread = new ThreadedPluginListener(theTP, event);
        thread.start();
    }
}
