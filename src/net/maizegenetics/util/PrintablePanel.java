/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */

//Title:        QTP Analyzer
//Version:
//Copyright:    Copyright (c) 1997
//Author:       Ed Buckler
//Company:      NCSU
//Description:  Your description

package net.maizegenetics.util;

//import com.sun.image.codec.jpeg.JPEGCodec;
//import com.sun.image.codec.jpeg.JPEGEncodeParam;
//import com.sun.image.codec.jpeg.JPEGImageEncoder;

import javax.imageio.ImageIO;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterJob;
import java.io.*;

public class PrintablePanel extends JPanel implements Printable {

  public PrintablePanel() {
    super();
  }


  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
  }

  public void jpegPrint(File saveFile) {
    if (saveFile== null) {return;}
    BufferedImage img;
    Dimension d=this.getSize();
    img = (BufferedImage)createImage(d.width, d.height);
    Graphics gbox=img.getGraphics();
    //This is where I used primitive graphic calls to draw the image.
    this.paintComponents(gbox);
    //Create a fileoutputstream.
    //FileOutputStream fos;
    try{
        
      // Deprecated code. should be removed.  
      //fos = new FileOutputStream(saveFile);
      //JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(fos);
      //JPEGEncodeParam params=encoder.getDefaultJPEGEncodeParam(img);
      //params.setQuality(1.0f,false);
      //encoder.encode(img,params);
      //fos.flush();
      //fos.close();
      
      ImageIO.write(img, "JPEG", saveFile);
    }
    catch(FileNotFoundException fe){System.out.println(fe);}
    catch(IOException ioe){System.out.println(ioe);}
  }

  public int print(Graphics g, PageFormat pf, int pageIndex) {
    if (pageIndex != 0) return NO_SUCH_PAGE;
    Graphics2D g2 = (Graphics2D)g;
    g2.translate(pf.getImageableX(), pf.getImageableY());
    double pageHeight = pf.getImageableHeight();
    double pageWidth = pf.getImageableWidth();
    double tableWidth = (double) this.getWidth();
    double tableHeight = (double) this.getHeight();
    double scaleW =  pageWidth / tableWidth;
    double scaleH =  pageHeight/tableHeight;
    double maxScale=(scaleW>scaleH)?scaleH:scaleW;
    System.out.println(scaleW+"=W  H="+scaleH+"   maxScale="+maxScale);
    g2.scale(maxScale,maxScale);
    this.paint(g2);
    return PAGE_EXISTS;
  }



  public void sendToPrinter() {
        PrinterJob printJob = PrinterJob.getPrinterJob();
        printJob.setPrintable(this);
        if (printJob.printDialog()) {
            try {
                printJob.print();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
  }


  public void svgSave(File saveFile) {
      if(saveFile==null) return;
      // Get a DOMImplementation
      DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();

      // Create an instance of org.w3c.dom.Document
      Document document = domImpl.createDocument(null, "svg", null);

      // Create an instance of the SVG Generator
      SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

      // Ask the test to render into the SVG Graphics2D implementation
      this.getComponent(0).paint(svgGenerator);

      // Finally, stream out SVG to the standard output using UTF-8
      // character to byte encoding
      boolean useCSS = true; // we want to use CSS style attribute
      try{
      FileOutputStream fos = new FileOutputStream(saveFile);
      Writer out = new OutputStreamWriter(fos, "UTF-8");
      svgGenerator.stream(out, useCSS);
      fos.flush();
      fos.close();
      }
      catch(Exception ee)
        {System.out.println("Error with svg button"+ee);}
  }

}