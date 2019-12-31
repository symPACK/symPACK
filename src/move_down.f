      subroutine  move_down
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
*
*       swap repeatedly with smallest among children
*       --------------------------------------------
*
        integer     fld_next1, fld_next2, fld_now, fld_next,
     &              val1, val2, val, vtx
*
        integer         infty
        parameter   (   infty = 2000000000 )
*       parameter   (   infty = 200000 )
*
        fld_now = fld_index
  100   continue
            if  ( fld_now .le. heap_size )  then
*
                val = heap(fld_now*2-1)
*
                fld_next1 = fld_now*2
                fld_next2 = fld_next1 + 1
*
                if  ( fld_next1 .le. heap_size )  then
                    val1 = heap(fld_next1*2-1)
                else
                    val1 = infty
                endif
*
                if  ( fld_next2 .le. heap_size )  then
                    val2 = heap(fld_next2*2-1)
                else
                    val2 = infty
                endif
*
                if  ( (val .le. val1) .and. (val .le. val2) )  then
                    fld_now = heap_size + 1
                else
                    if  ( val1 .le. val2 )  then
                        fld_next = fld_next1
                    else
                        fld_next = fld_next2
                    endif
*
*                   ----
*                   swap
*                   ----
*
                    fld_next1 = fld_now*2
                    fld_next2 = fld_next*2
*
                    vtx = heap(fld_next1)
*
                    vtx_to_heap(vtx) = fld_next
                    vtx_to_heap(heap(fld_next2)) = fld_now
*
                    heap(fld_next1) = heap(fld_next2)
                    heap(fld_next2) = vtx
*
                    fld_next1 = fld_next1 - 1
                    fld_next2 = fld_next2 - 1
                    heap(fld_next1) = heap(fld_next2)
                    heap(fld_next2) = val
                    fld_now = fld_next
                endif
                go to 100
*
            endif
        return
*
      end
