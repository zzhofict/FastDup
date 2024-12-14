/* yarn.h -- generic interface for thread operations
 * Copyright (C) 2008, 2011, 2012, 2015, 2018, 2019, 2020 Mark Adler
 * Version 1.7  12 Apr 2020  Mark Adler
 */

/*
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the author be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Mark Adler
  madler@alumni.caltech.edu
 */

/* Basic thread operations

   This interface isolates the local operating system implementation of threads
   from the application in order to facilitate platform independent use of
   threads.  All of the implementation details are deliberately hidden.

   Assuming adequate system resources and proper use, none of these functions
   can fail.  As a result, any errors encountered will cause an exit() to be
   executed, or the execution of your own optionally-provided abort function.

   These functions allow the simple launching and joining of threads, and the
   locking of objects and synchronization of changes of objects.  The latter is
   implemented with a single lock_t type that contains an integer value.  The
   value can be ignored for simple exclusive access to an object, or the value
   can be used to signal and wait for changes to an object.

   -- Arguments --

   thread *thread;          identifier for launched thread, used by join
   void probe(void *);      pointer to function "probe", run when thread starts
   void *payload;           single argument passed to the probe function
   lock_t *lock_t;              a lock_t with a value -- used for exclusive access to
                            an object and to synchronize threads waiting for
                            changes to an object
   long val;                value to set lock_t, increment lock_t, or wait for
   int n;                   number of threads joined

   -- Thread functions --

   thread = launch(probe, payload) - launch a thread -- exit via probe() return
   join(thread) - join a thread and by joining end it, waiting for the thread
        to exit if it hasn't already -- will free the resources allocated by
        launch() (don't try to join the same thread more than once)
   n = join_all() - join all threads launched by launch() that are not joined
        yet and free the resources allocated by the launches, usually to clean
        up when the thread processing is done -- join_all() returns an int with
        the count of the number of threads joined (join_all() should only be
        called from the main thread, and should only be called after any calls
        of join() have completed)

   -- Lock functions --

   lock_t = new_lock(val) - create a new lock_t with initial value val (lock_t is
        created in the released state)
   possess(lock_t) - acquire exclusive possession of a lock_t, waiting if necessary
   twist(lock_t, [TO | BY], val) - set lock_t to or increment lock_t by val, signal
        all threads waiting on this lock_t and then release the lock_t -- must
        possess the lock_t before calling (twist releases, so don't do a
        release() after a twist() on the same lock_t)
   wait_for(lock_t, [TO_BE | NOT_TO_BE | TO_BE_MORE_THAN | TO_BE_LESS_THAN], val)
        - wait on lock_t value to be, not to be, be greater than, or be less than
        val -- must possess the lock_t before calling, will possess the lock_t on
        return but the lock_t is released while waiting to permit other threads
        to use twist() to change the value and signal the change (so make sure
        that the object is in a usable state when waiting)
   release(lock_t) - release a possessed lock_t (do not try to release a lock_t that
        the current thread does not possess)
   val = peek_lock(lock_t) - return the value of the lock_t (assumes that lock_t is
        already possessed, no possess or release is done by peek_lock())
   free_lock(lock_t) - free the resources allocated by new_lock() (application
        must assure that the lock_t is released before calling free_lock())

   -- Memory allocation ---

   yarn_mem(better_malloc, better_free) - set the memory allocation and free
        routines for use by the yarn routines where the supplied routines have
        the same interface and operation as malloc() and free(), and may be
        provided in order to supply thread-safe memory allocation routines or
        for any other reason -- by default malloc() and free() will be used

   -- Error control --

   yarn_prefix - a char pointer to a string that will be the prefix for any
        error messages that these routines generate before exiting -- if not
        changed by the application, "yarn" will be used
   yarn_abort - an external function that will be executed when there is an
        internal yarn error, due to out of memory or misuse -- this function
        may exit to abort the application, or if it returns, the yarn error
        handler will exit (set to NULL by default for no action)
 */
#pragma once

#include <stdlib.h>

namespace yarn {
    
extern char *yarn_prefix;
extern void (*yarn_abort)(int);

void yarn_mem(void *(*)(size_t), void (*)(void *));

typedef struct thread_s thread;
thread *launch_(void (*)(void *), void *, char const *, long);
#define LAUNCH(a, b) launch_(a, b, __FILE__, __LINE__)
void join_(thread *, char const *, long);
#define JOIN(a) join_(a, __FILE__, __LINE__)
int join_all_(char const *, long);
#define JOIN_ALL() join_all_(__FILE__, __LINE__)

typedef struct lock_s lock_t;
lock_t *new_lock_(long, char const *, long);
#define NEW_LOCK(a) new_lock_(a, __FILE__, __LINE__)
void possess_(lock_t *, char const *, long);
#define POSSESS(a) possess_(a, __FILE__, __LINE__)
void release_(lock_t *, char const *, long);
// #define release(a) release_(a, __FILE__, __LINE__)
#define RELEASE(a) release_(a, __FILE__, __LINE__)
enum twist_op { TO, BY };
void twist_(lock_t *, enum twist_op, long, char const *, long);
#define TWIST(a, b, c) twist_(a, b, c, __FILE__, __LINE__)
enum wait_op {
    TO_BE,
    /* or */ NOT_TO_BE, /* that is the question */
    TO_BE_MORE_THAN,
    TO_BE_LESS_THAN
};
void wait_for_(lock_t *, enum wait_op, long, char const *, long);
#define WAIT_FOR(a, b, c) wait_for_(a, b, c, __FILE__, __LINE__)
long peek_lock(lock_t *);
#define PEEK_LOCK(a) peek_lock(a)
void free_lock_(lock_t *, char const *, long);
#define FREE_LOCK(a) free_lock_(a, __FILE__, __LINE__)

};  // namespace yarn
