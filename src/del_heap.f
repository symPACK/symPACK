      subroutine  del_heap
     &      ( heap, heap_size, vtx_to_heap, vtx )
*
        integer     heap(*)
        integer     heap_size
        integer     vtx_to_heap(*)
        integer     vtx
*
        integer     i, j, k, v, val
*
        v = vtx
        if  ( heap_size .eq. 0 )  return
        i = vtx_to_heap(v)
        if  ( i .lt. heap_size )  then
            vtx_to_heap(v) = 0
            j = heap_size*2
            k = i*2
            v = heap(j)
            heap(k) = v
            val = heap(j-1)
            vtx_to_heap(v) = i
            heap_size = heap_size - 1
            call  mod_heap ( heap, heap_size, vtx_to_heap, v, val )
        else
            vtx_to_heap(v) = 0
            heap_size = heap_size - 1
        endif
        return
*
      end
