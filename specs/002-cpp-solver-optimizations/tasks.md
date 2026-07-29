# Tasks: C++ Solver Performance Optimizations

**Input**: Design documents from `specs/002-cpp-solver-optimizations/`

**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/library_api.md, quickstart.md

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

---

## Phase 1: Setup

**Purpose**: Project preparation for optimizations. (No setup tasks required beyond existing headers).

---

## Phase 2: Foundational

**Purpose**: Define the persistent workspace structures required by all downstream tasks.

- [x] T001 [TDD] [SUBAGENT] Define `CloudJ::RadiativeSolver::Workspace` struct with flat vectors and `resize` method in `include/cloudj/radiative_solver.hpp`

---

## Phase 3: User Story 1 - Zero-Allocation Column Solves (Priority: P1) 🎯 MVP

**Goal**: Eliminate dynamic memory allocation inside the calculation loop.

- [x] T002 [P] [US1] Update `CloudJ::Engine` class to instantiate and hold a persistent `RadiativeSolver::Workspace` member in `include/cloudj/cloudj.hpp`
- [x] T003 [P] [US1] Modify `CloudJ::Engine::calculate_photolysis_rates` to call `solver_ws.resize(lu)` prior to calculations in `include/cloudj/cloudj.hpp`
- [ ] T004 [US1] Update `RadiativeSolver::BLKSLV` method signature to accept `Workspace& ws` parameter in `include/cloudj/radiative_solver.hpp`
- [ ] T005 [US1] Refactor `BLKSLV` to use the provided `ws` arrays instead of allocating local `std::vector` buffers in `include/cloudj/radiative_solver.hpp`
- [ ] T006 [REVIEW] [US1] Pass the engine's `solver_ws` down into `BLKSLV` invocations from `include/cloudj/cloudj.hpp` (Note: `BLKSLV` is called indirectly, trace the call graph to ensure it's passed correctly, likely through `Photolysis::JRATET` or similar if it's nested).

**Checkpoint**: At this point, dynamic heap allocations during 10,000 sequential column solves should be exactly zero after the engine is initialized.

---

## Phase 4: User Story 2 - SIMD Compiler Vectorization (Priority: P2)

**Goal**: Unroll 4x4 matrix operations for SIMD efficiency.

- [x] T007 [P] [SUBAGENT] [US2] Unroll initialization loops for `b`, `aa`, `cc`, `a`, `c`, `h` matrices inside `GEN_ID` in `include/cloudj/radiative_solver.hpp`
- [x] T008 [P] [SUBAGENT] [US2] Unroll 4x4 nested matrix multiplication and summation loops inside `GEN_ID` in `include/cloudj/radiative_solver.hpp`
- [x] T009 [P] [SUBAGENT] [US2] Ensure `solve_lu_4x4` is fully unrolled (already manually unrolled in previous implementation, verify and optimize if needed) in `include/cloudj/radiative_solver.hpp`
- [x] T010 [P] [SUBAGENT] [US2] Unroll back-substitution and mean J/H calculation loops inside `BLKSLV` in `include/cloudj/radiative_solver.hpp`

**Checkpoint**: At this point, loop boundary checking branches for M=4 should be eliminated.

---

## Phase 5: User Story 3 - Contiguous Cache Line Hits for Actinic Flux (Priority: P1)

**Goal**: Flatten the actinic flux `fff` array for better cache locality.

- [x] T011 [P] [US3] Change `fff` from `std::vector<std::vector<double>>` to a flat `std::vector<double>` in `include/cloudj/cloudj.hpp`
- [x] T012 [P] [US3] Create an `mdspan_2d_mut` view (`layout_left`) over the flat `fff` array in `include/cloudj/cloudj.hpp`
- [x] T013 [P] [US3] Update `Photolysis::JRATET` signature to accept `mdspan_2d_mut fff` instead of `std::vector<std::vector<double>>&` in `include/cloudj/photolysis.hpp`
- [x] T014 [US3] Refactor `Photolysis::JRATET` to use `mdspan` 2D indexing (`fff(k, l)`) instead of nested vector indexing (`fff[k][l]`) in `include/cloudj/photolysis.hpp`

**Checkpoint**: Actinic flux memory accesses sequentially scan adjacent indices.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Validation and benchmarking.

- [x] T015 [TDD] Compile and run unit tests (`test_library_api.cpp`) to verify API changes and numerical parity.
- [x] T016 Run the 500,000-run property fuzzer to guarantee $10^{-12}$ exact numerical parity remains intact.
- [x] T017 Run the high-resolution benchmarking script (`run_benchmarks.py`) to verify the $\ge 1.1x$ in-memory throughput speedup.

---

## Dependencies & Execution Order

1. **Phase 1 (Setup)**: N/A
2. **Phase 2 (Foundational)**: Prerequisite for Phase 3.
3. **Phase 3 (User Story 1)**: Modifies core solver signatures, should be done first to establish the workspace pattern.
4. **Phase 4 (User Story 2)**: Can be done independently of US1 and US3, but modifies the same `radiative_solver.hpp` file, so execute sequentially after Phase 3.
5. **Phase 5 (User Story 3)**: Modifies `cloudj.hpp` and `photolysis.hpp`. Can be executed in parallel with Phase 4 if file conflicts are managed, but sequential is safer.
6. **Phase 6 (Polish)**: Depends on all prior phases.

## Implementation Strategy

### Incremental Delivery
1. Implement persistent workspaces (Phase 2 & 3) and verify zero allocations.
2. Flatten actinic flux array (Phase 5) and verify cache-locality speedup.
3. Unroll 4x4 loops (Phase 4) and verify SIMD speedup.
4. Final extensive validation (Phase 6) to guarantee no loss of numerical parity.
