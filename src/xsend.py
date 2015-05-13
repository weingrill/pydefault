#!/usr/bin/python
# -*- coding: utf-8 -*-
# $Id: xsend.py,v 1.8 2006/10/06 12:30:42 normanr Exp $
'''
Created on May 7, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
#!/usr/bin/python
import sys, xmpp, time

tojid='jweingrill@communigate.aip.de'
text='This is a Jabber test'

jidparams={'jid':'jweingrill@communigate.aip.de', 'password':'jo2712we'}

jid = xmpp.protocol.JID(jidparams['jid'])
cl = xmpp.Client(jid.getDomain(),debug=[])

con = cl.connect()
if not con:
    print 'could not connect!'
    sys.exit()
    
print 'connected with',con
auth = cl.auth(jid.getNode(),jidparams['password'],resource=jid.getResource())

if not auth:
    print 'could not authenticate!'
    sys.exit()
    
print 'authenticated using',auth

#cl.SendInitPresence(requestRoster=0)   # you may need to uncomment this for old server
ident = cl.send(xmpp.protocol.Message(tojid,text))
print 'sent message with id', ident

time.sleep(1)   # some older servers will not send the message if you disconnect immediately after sending

#cl.disconnect()