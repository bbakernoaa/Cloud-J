# Tasks: Advanced C++ Solver Performance Optimizations

**Input**: Design documents from `specs/003-advanced-performance-optimizations/`

**Prerequisites**: plan.md (required), spec.md (required for user stories)

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

---

## Phase 1: Setup

**Purpose**: Project preparation for optimizations. (No setup tasks required beyond existing headers).

---

## Phase 2: Foundational

**Purpose**: Apply struct alignment to prevent False Sharing.

- [x] T001 [P] [US2] Apply C++11 standard alignment specifyer `alignas(64)` to the `Workspace` struct declaration in `include/cloudj/radiative_solver.hpp`

---

## Phase 3: User Story 1 - Division-Free Cross-Section Interpolation (Priority: P1) 🎯 MVP

**Goal**: Eliminate FDIV operations from the inner interpolation loop.

- [x] T002 [US1] Pre-compute reciprocal temperature span intervals (`inv_t12` and `inv_t23`) during species loading inside `SpecData` inside `include/cloudj/photolysis.hpp`
- [x] T003 [US1] Modify `CrossSections::interpolate` in `include/cloudj/cross_sections.hpp` to accept pre-computed reciprocal variables and replace divisions with FMA multiplications.

---

## Phase 4: User Story 3 - Algebraic Reciprocal manual 4x4 LU Solves (Priority: P1)

**Goal**: Replace divisions inside the manual 4x4 LU solver with multiplications.

- [x] T004 [US3] Refactor `solve_lu_4x4` in `include/cloudj/radiative_solver.hpp` to pre-calculate pivot reciprocals and multiply instead of divide for all dependent elements.

---

## Phase 5: Polish & Cross-Cutting Concerns

**Purpose**: Validation and benchmarking.

- [x] T005 [TDD] Compile and run unit tests (`test_library_api.cpp`) to verify API changes and numerical parity.
- [x] T006 Run the 500,000-run property fuzzer to guarantee $10^{-12}$ exact numerical parity remains intact.
- [x] T007 Run the high-resolution benchmarking script (`run_benchmarks.py`) to verify the speedups.

---

## Dependencies & Execution Order

1. **Phase 1 (Setup)**: N/A
2. **Phase 2 (Foundational)**: Struct alignment is completely independent.
3. **Phase 3 (User Story 1)**: Modifies photolysis and cross-section parameters.
4. **Phase 4 (User Story 3)**: Modifies the 4x4 LU solver.
5. **Phase 5 (Polish)**: Depends on all prior phases.
