CONTIKI = ../../contiki-3.0

all: link_test

CONTIKI_WITH_RIME = 1
DEFINES+=NETSTACK_CONF_RDC=nullrdc_driver
DEFINES+=NETSTACK_CONF_MAC=nullmac_driver
DEFINES+=NETSTACK_CONF_FRAMER=framer_802154
include $(CONTIKI)/Makefile.include
