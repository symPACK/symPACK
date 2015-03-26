/**
 * @file test_list.c
 * @author Mark Hoemmen
 * @since 11 July 2006
 * @date Time-stamp: <2008-07-16 10:21:46 mhoemmen>
 * 
 * Correctness tests for the list_t data structure.  Executable
 * returns 0 if all the tests pass, otherwise it breaks with an
 * assertion.
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_includes.h"
#include <bebop/util/list.h>

int* 
new_int (const int i)
{
  int* p = bebop_malloc (sizeof (int));
  *p = i;
  return p;
}

int
main (int argc, char** argv)
{
  int* item = NULL;
  int* item2 = NULL;
  int* item3 = NULL;
  list_t L;
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);
  bebop_log (1, "test_list: starting tests\n");
  
  L = list_create ();
  assert (list_empty_p (L));
  assert (L.head == NULL);
  assert (L.tail == NULL);

  item = new_int (42);
  assert (*item == 42);
  L = list_append_item (L, item);
  assert (! list_empty_p (L));
  assert (L.head != NULL);
  assert ((L.head)->item == item);
  assert (list_car (L) == item);
  assert (L.head == L.tail);
  assert ((L.head)->item != NULL);
  assert (list_car (L) != NULL);
  assert (*( (int*) ((L.head)->item) ) == 42);
  assert ( *( (int*) list_car (L) ) == 42 );
  assert ((L.tail)->next == NULL);
  assert ((list_cdr (L)).head == NULL);
  assert ((list_cdr (L)).tail == NULL);

  item2 = new_int (43);
  assert (*item2 == 43);
  L = list_append_item (L, item2);
  assert (! list_empty_p (L));
  assert ((L.head)->item == item);
  assert ((L.tail)->item == item2);
  assert (L.head != L.tail);
  assert (*( (int*) ((L.head)->item) ) == 42);
  assert (*( (int*) ((L.tail)->item) ) == 43);
  assert ((L.head)->next == L.tail);
  assert ((L.tail)->next == NULL);

  item3 = new_int (44);
  assert (*item3 == 44);
  L = list_append_item (L, item3);
  assert (! list_empty_p (L));
  assert ((L.head)->item == item);
  assert ((L.head)->next->item == item2);
  assert ((L.head)->next->next->item == item3);
  assert ((L.head)->next->next == L.tail);
  assert (L.head != L.tail);
  assert ((L.head)->next != L.tail);
  assert (*( (int*) ((L.head)->item) ) == 42);
  assert (*( (int*) ((L.tail)->item) ) == 44);
  assert ((L.head)->next->next->next == NULL);
  
  L = list_destroy (L);
  assert (L.head == NULL);
  assert (L.tail == NULL);

  bebop_log (1, "test_list: passed all tests\n");
  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS; /* to pacify the compiler */
}

