##############################################################################

#MK_TOP = ../../../..
MK_TOP  = /export/ciao_from_source/ciao-4.8/src/
KJG = /export/ciao


include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmradar
LIB_FILES = 
PAR_FILES         = dmradar.par
INC_FILES         = 
XML_FILES         = dmnautilus.xml

SRCS	=           dmnautilus.c t_dmnautilus.c

LOCAL_LIBS = -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I$(MK_TOP)/da/analysis/dmtools/dmimgio/

OBJS = $(SRCS:.c=.o)


MAKETEST_SCRIPT   = dmnautilus.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------



$(EXEC): $(OBJS)
	$(LINK)
	@echo


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                Building dmnautilus                       | "
	@echo "   \----------------------------------------------------------/ "


kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)

