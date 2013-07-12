package net.maizegenetics.util;

import java.util.*;
import javax.mail.*;
import javax.mail.internet.*; 
import javax.activation.*;


/**
 * This class can be used to send email.  When the host is "appsmtp.mail.cornell.edu", any machine within the Cornell
 * network can be used to send an email.
 * Reference: https://confluence.cornell.edu/login.action?os_destination=https%3A%2F%2Fconfluence.cornell.edu%2Fdisplay%2Fcitapps%2FHow%2Bto%2Bconfigure%2Bapplications%2Bto%2Bsend%2Bmail%2Bthrough%2Bthe%2BCornell%2Bmail%2Bservers.
 */
public class SMTPClient { 
    
    private String host;
    private String fromAddress;
    private String toAddress;
    private static MimeMessage message; 
    
    public SMTPClient(String host, String toAddress){
        
        Properties properties = System.getProperties();
        properties.setProperty("mail.smtp.host", host);
        
        // Get the default Session object.
        Session session = Session.getDefaultInstance(properties);
        message = new MimeMessage(session); 

        try{
            message.addRecipient(Message.RecipientType.TO, new InternetAddress(toAddress));
            message.setFrom(new InternetAddress(toAddress));
        }catch(javax.mail.internet.AddressException ae){ /* ignore */  }
         catch(javax.mail.MessagingException me){  /* ignore */  }
    }

    public void sendMessageWithAttachment(String subject, String msg, String fileAttachment) throws javax.mail.MessagingException{
        
        MimeBodyPart messageBodyPart = new MimeBodyPart();
        messageBodyPart.setText(msg);
        
        Multipart multipart = new MimeMultipart();
        multipart.addBodyPart(messageBodyPart);

        // Part two is attachment
        messageBodyPart = new MimeBodyPart();
        DataSource source = new FileDataSource(fileAttachment);
        messageBodyPart.setDataHandler(new DataHandler(source));
        messageBodyPart.setFileName(fileAttachment);
        multipart.addBodyPart(messageBodyPart);

        // Put parts in message
        message.setContent(multipart);

        Transport.send(message);
        
    }

    public void sendMessage(String subject, String msg) throws javax.mail.MessagingException{
        message.setSubject(subject);
        message.setText(msg);
        Transport.send(message);
    }
}
