void
add_partlineages426 (long numpop, timelist_fmt ** timevector, world_fmt* world)
{
    // this points to the MRCA
    long T;
    long i, pop;
    vtlist *tl = (*timevector)->tl;
    vtlist *tli;
    vtlist *tli1;
    long from;
    long to;
    long *lineages;
    char type;
    T = (*timevector)->T;
    T = T - 2;
    memset(tl[T+1].lineages, 0, (size_t) numpop * sizeof(long));
    memset(tl[T].lineages, 0, (size_t) numpop * sizeof(long));
    tl[T+1].lineages[tl[T].eventnode->pop] = 1;
    //tl[T].lineages[tl[T].eventnode->pop] = 2;
    for (i = T; i >= 0; i--)
    {
        tli1 = &tl[i+1];
        tli = &tl[i];
        lineages = tli->lineages;
        for(pop=0;pop<numpop;pop++)
            lineages[pop] = tli1->lineages[pop];
        //memcpy(tli->lineages,tli1->lineages, (size_t) numpop * sizeof(long));
        from = tli->from;
        to = tli->to;
        //if(tli->eventnode != NULL)
        type = tli->eventnode->type;
        //else
        //type = 'b'; // boundary
        
        switch(type)
        {
            case 't':
                lineages[to] -= 1;
                if (lineages[to] < 0)
                {
                    print_dead_tree(numpop, *timevector, world);
                    error("lineages < 0 in addpartlineages");
                }
                break;
            case 'i':
                lineages[to] += 1;
                break;
            case 'm':
            case 'd':
                lineages[to] += 1;
                lineages[from] -= 1;
                if (lineages[from] < 0)
                {
                    print_dead_tree(numpop, *timevector, world);
                    error("lineages < 0 in addpartlineages");
                }
                break;
            case 'r':
                //case 'b':
                warning("Root node in tymelist, this should not happen");
            default:
            {
                printf("%i> failed in constructing timelist: with type=%i\n",myID,type);
                print_dead_tree(numpop, *timevector, world);
                error("funny node received and died in add_partlineages");
            }
                //  break;
        }
    }
}

