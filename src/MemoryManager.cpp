/**
 * A global memory manager that allows allocation of
 * memory blocks that are shared across threads and
 * processes. Memory blocks are identified by strings.
 */

#include <map>
#include <omp.h>
#include <string>
#include "config.h"
#include "MemoryManager.h"

using namespace std;
using namespace MemoryManager;

namespace MemoryManager {
    map<string, struct memory_block*> blocks;

    /**
     * Allocate a memory block with the given name
     * (if the block has not already been allocated).
     *
     * name: Name of block to allocate memory for.
     */
    void *allocate(const string &name, const size_t size) {
        #pragma omp critical (MemoryManager_allocate)
        {
            if (!block_exists(name)) {
                blocks[name] = new struct memory_block;

                blocks[name]->name = name;
                blocks[name]->size = size;
                blocks[name]->ptr  = malloc(size);
            }
        }

        return blocks[name]->ptr;
    }

    /**
     * Check if a memory block with the given name exists.
     * 
     * name: Name of block to check for.
     */
    bool block_exists(const string &name) {
        return (blocks.find(name) != blocks.end());
    }

    /**
     * Deallocate memory block.
     *
     * name: Name of block to deallocate.
     */
    void deallocate(const string &name) {
        free(blocks[name]->ptr);
        delete blocks[name];

        blocks.erase(name);
    }
}
