      subroutine  ins_heap
     &      ( heap, heap_size, vtx_to_heap, vtx, val )
*
        integer     heap(*)
        integer     heap_size
        integer     vtx_to_heap(*)
        integer     vtx
        integer     val
*
        integer     i, j
*
        heap_size = heap_size + 1
        i = heap_size
        vtx_to_heap(vtx) = i
        j = i*2
        heap(j) = vtx
        heap(j-1) = val
        if  ( i .eq. 1 )  return
        call  move_up ( heap, heap_size, i, vtx_to_heap )
        return
*
      end
