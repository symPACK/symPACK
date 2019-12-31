      subroutine  move_up
     &      ( heap, heap_size, fld_index, vtx_to_heap )
*
        integer     heap(*)
        integer     heap_size
        integer     fld_index
        integer     vtx_to_heap(*)
*
*       --------------------------------------------
*       field at 2*fld_index-1, 2*fld_index
*                   value          vtx
*
*       children at 2*fld_index, 2*fld_index +1
*       parent   at fld_index/2 
*
*       swap repeatedly with with parnet if value at
*       at parent is larger
*       --------------------------------------------
*
        integer     fld_next1, fld_next2, fld_now, fld_next,
     &              val1, val, vtx
*
        fld_now = fld_index
  100   continue
            if  ( fld_now .gt. 1 )  then
*
                fld_next = fld_now/2
*
                fld_next1 = fld_now*2
                val = heap(fld_next1-1)
*
                fld_next2 = fld_next*2
                val1 = heap(fld_next2-1)
*
                if  ( val .ge. val1 )  then
                    fld_now = 0
                else
*
*                   ----------
*                   val < val1
*
*                   swap ...
*                   ----------
*
                    vtx = heap(fld_next1)
                    vtx_to_heap(vtx) = fld_next
                    vtx_to_heap(heap(fld_next2)) = fld_now
*
                    heap(fld_next1) = heap(fld_next2)
                    heap(fld_next2) = vtx
*
                    fld_next1 = fld_next1 - 1
                    fld_next2 = fld_next2 - 1
                    heap(fld_next1) = val1
                    heap(fld_next2) = val
                    fld_now = fld_next
                endif
                go to 100
            endif
        return
*
      end
