#ifndef _MEMORY_MANAGER_H
#define _MEMORY_MANAGER_H

#include <string>

namespace MemoryManager {
    struct memory_block {
        std::string name;
        size_t size;
        void *ptr;
    };
    
    void *allocate(const std::string&, const size_t);
    bool block_exists(const std::string&);
    void deallocate(const std::string&);
    void *get_block(const std::string&);
    void *set_block(const std::string&, const size_t, void*);
}

#endif/*_MEMORY_MANAGER_H*/
