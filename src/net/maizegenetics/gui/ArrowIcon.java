/*
 * ArrowIcon.java
 *
 * Created on February 4, 2009, 5:57 PM
 *
 */
package net.maizegenetics.gui;

import javax.swing.Icon;
import java.awt.*;

public class ArrowIcon implements Icon {

    public static final int UP = 0;
    public static final int DOWN = 1;
    public static final int LEFT = 2;
    public static final int RIGHT = 3;
    private int direction;
    private Polygon pagePolygon = new Polygon(new int[]{2, 4, 4, 10, 10, 2},
            new int[]{4, 4, 2, 2, 12, 12},
            6);
    private int[] arrowX = {4, 9, 6};
    private Polygon arrowUpPolygon =
            new Polygon(arrowX, new int[]{10, 10, 4}, 3);
    private Polygon arrowDownPolygon =
            new Polygon(arrowX, new int[]{6, 6, 11}, 3);
    private int[] arrowY = {5, 10, 7};
    private Polygon arrowLeftPolygon =
            new Polygon(new int[]{9, 9, 3}, arrowY, 3);
    private Polygon arrowRightPolygon =
            new Polygon(new int[]{4, 4, 9}, arrowY, 3);

    public ArrowIcon(int which) {
        direction = which;
    }

    public int getIconWidth() {
        return 14;
    }

    public int getIconHeight() {
        return 14;
    }

    public void paintIcon(Component c, Graphics g, int x, int y) {

        g.setColor(Color.black);
        pagePolygon.translate(x, y);
        g.drawPolygon(pagePolygon);
        pagePolygon.translate(-x, -y);
        if (direction == UP) {
            arrowUpPolygon.translate(x, y);
            g.fillPolygon(arrowUpPolygon);
            arrowUpPolygon.translate(-x, -y);
        } else if (direction == DOWN) {
            arrowDownPolygon.translate(x, y);
            g.fillPolygon(arrowDownPolygon);
            arrowDownPolygon.translate(-x, -y);
        } else if (direction == LEFT) {
            arrowLeftPolygon.translate(x, y);
            g.fillPolygon(arrowLeftPolygon);
            arrowLeftPolygon.translate(-x, -y);
        } else if (direction == RIGHT) {
            arrowRightPolygon.translate(x, y);
            g.fillPolygon(arrowRightPolygon);
            arrowRightPolygon.translate(-x, -y);
        } else {
            throw new IllegalArgumentException("ArrowIcon: paintIcon: unknown direction.");
        }

    }
}

