//
// Link check of the run side: this program is linked from every object of opalx_run and its
// third-party libraries, without opalx. A run-side object that refers to an input-side or bridge
// symbol makes the link fail with an undefined symbol; running the program only has to succeed.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
class Inform;

// Run-side files write to gmsg, which changes.cpp defines. Reading it makes a shared opalx_run a
// needed library of this program, so the linker also checks the undefined symbols of that library.
extern Inform* gmsg;

int main() {
    Inform* volatile message = gmsg;
    static_cast<void>(message);
    return 0;
}
