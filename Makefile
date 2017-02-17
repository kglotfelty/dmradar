##############################################################################

#MK_TOP = ../../../..
MK_TOP  = /export/ciao_from_source/ciao-4.9/src/
KJG = /export/ciao


include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmradar
LIB_FILES = 
PAR_FILES         = dmradar.par
INC_FILES         = 
XML_FILES         = dmradar.xml

SRCS	=           dmradar.c t_dmradar.c

LOCAL_LIBS = -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I$(MK_TOP)/da/analysis/dmtools/dmimgio/

OBJS = $(SRCS:.c=.o)


MAKETEST_SCRIPT   = dmradar.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------



$(EXEC): $(OBJS)
	$(LINK)
	@echo


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                Building dmradar                       | "
	@echo "   \----------------------------------------------------------/ "


kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)

