/*
 * Copyright (c) 2007, Swedish Institute of Computer Science.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the Institute nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE INSTITUTE AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE INSTITUTE OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *
 */

/**
 * \file
 *         Based on "Best-effort single-hop unicast example" from CONTIKI
 *         The goal is to measure link quality variation between two nodes
 *         Should use nullrdc_driver, nullmac_driver and framer_802154
 *         Button click to start
  */

#include "contiki.h"
#include "net/rime/rime.h"

#include "dev/cc2420/cc2420.h"
#include "dev/button-sensor.h"
#include "dev/leds.h"

#include <stdio.h>
#include <string.h>

// Configurations
#define RSSI_OFFSET 100
#define TRANS_INTER (CLOCK_SECOND/10)

// Global variables
static linkaddr_t dest_addr = {{0, 0}};
static uint8_t initialized = 0;
static uint8_t clicked = 0;
static process_event_t myev_num;
static int counter = 0;

/*---------------------------------------------------------------------------*/
PROCESS(link_test_process, "Link test");
AUTOSTART_PROCESSES(&link_test_process);
/*---------------------------------------------------------------------------*/
static void
recv_uc(struct unicast_conn *c, const linkaddr_t *from)
{
  printf("%X.%X>%X.%X|%sL%dR%d",
  from->u8[0], from->u8[1],
   linkaddr_node_addr.u8[0], linkaddr_node_addr.u8[1],
	 (char*) packetbuf_dataptr(),
   packetbuf_attr(PACKETBUF_ATTR_LINK_QUALITY),
   packetbuf_attr(PACKETBUF_ATTR_RSSI)+RSSI_OFFSET);
  printf("\n");

  if (!initialized) {
    printf("INIT %d\n", counter);
    initialized = 1;
    linkaddr_copy(&dest_addr, from);
  }
  process_post(&link_test_process, myev_num, 0);
}
/*---------------------------------------------------------------------------*/
static const struct broadcast_callbacks broadcast_call = \
  {(void (*)(struct broadcast_conn *, const linkaddr_t *)) recv_uc};
static const struct unicast_callbacks unicast_callbacks = {recv_uc};
static struct broadcast_conn broadcast;
static struct unicast_conn uc;
/*---------------------------------------------------------------------------*/
PROCESS_THREAD(link_test_process, ev, data)
{
  PROCESS_EXITHANDLER(unicast_close(&uc);)

  PROCESS_BEGIN();
  static struct etimer et;

  SENSORS_ACTIVATE(button_sensor);
  // To change radio TX power
  // cc2420_set_txpower(3);

  myev_num = process_alloc_event();
  unicast_open(&uc, 146, &unicast_callbacks);
  broadcast_open(&broadcast, 129, &broadcast_call);

  printf("My addr: %X.%X\n",linkaddr_node_addr.u8[0],linkaddr_node_addr.u8[1]);

  while(1) {
    char message[10];

    if(initialized || clicked) {
      etimer_set(&et, TRANS_INTER);
    }

    PROCESS_WAIT_EVENT();

    if (ev == sensors_event && data == &button_sensor) {
      printf("BTN\n");
      clicked = !clicked;
    }

    if (ev == myev_num) {
      etimer_set(&et, TRANS_INTER/2);
      PROCESS_WAIT_EVENT();
    }
    if (etimer_expired(&et)) {
      sprintf(message, "%4d", ++counter);

      packetbuf_copyfrom(message, strlen(message)+1);
      // To change packet payload size
      // packetbuf_set_datalen(50);

      // printf("%X.%X>%X.%X|%s",
      //   linkaddr_node_addr.u8[0], linkaddr_node_addr.u8[1],
      //   dest_addr.u8[0], dest_addr.u8[1],
      //   (char*) packetbuf_dataptr());
      //   printf("\n");
      if (linkaddr_cmp(&dest_addr, &linkaddr_null)) {
        broadcast_send(&broadcast);
      } else {
        unicast_send(&uc, &dest_addr);
      }
    }

  }

  PROCESS_END();
}
/*---------------------------------------------------------------------------*/
