/*
 * AlignmentTableCellRenderer
 */
package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Component;

import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

import net.maizegenetics.pal.alignment.AlignmentMask;

/**
 *
 * @author terry
 */
public class AlignmentTableCellRenderer extends DefaultTableCellRenderer {

    private final AlignmentTableModel myAlignmentTableModel;
    private final AlignmentMask[] myMask;

    public AlignmentTableCellRenderer(AlignmentTableModel model) {
        myAlignmentTableModel = model;
        myMask = null;
    }

    public AlignmentTableCellRenderer(AlignmentTableModel model, AlignmentMask mask) {
        this(model, new AlignmentMask[]{mask});
    }

    public AlignmentTableCellRenderer(AlignmentTableModel model, AlignmentMask[] mask) {
        myAlignmentTableModel = model;
        myMask = mask;
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        Color theColor = null;
        if ((myMask == null) || (myMask.length == 0)) {
            return comp;
        } else if (myMask.length == 1) {

            if (myMask[0].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                theColor = myMask[0].getColor();
            }

        } else {

            int red = 0;
            int green = 0;
            int blue = 0;

            boolean changed = false;
            for (int i = 0; i < myMask.length; i++) {
                if (myMask[i].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                    red = red + myMask[i].getColor().getRed();
                    green = green + myMask[i].getColor().getGreen();
                    blue = blue + myMask[i].getColor().getBlue();
                    changed = true;
                }
            }

            if (changed) {
                red = red % 256;
                green = green % 256;
                blue = blue % 256;
                theColor = new Color(red, green, blue);
            }

        }

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else {
            comp.setBackground(theColor);
        }

        return comp;

    }
}
