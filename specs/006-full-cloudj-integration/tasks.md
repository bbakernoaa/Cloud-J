# Tasks: Full Cloud-J Integration and Complete Porting

**Input**: Design documents from `specs/006-full-cloudj-integration/`

**Prerequisites**: plan.md (required), spec.md (required for user stories)

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

---

## Phase 1: Setup

**Purpose**: N/A (No setup tasks required beyond existing headers).

---

## Phase 2: Foundational

**Purpose**: Port Legendre expansions and scattered flux orchestrators.

- [x] T001 Port the `LEGND0` ordinary Legendre polynomial expansion subroutine from `cldj_fjx_sub_mod.F90` to `include/cloudj/radiative_solver.hpp`.
- [x] T002 Port the `MIESCT` scattered flux orchestrator from `cldj_fjx_sub_mod.F90` to `include/cloudj/radiative_solver.hpp`, wiring up Legendre expansions and block solver calls (`BLKSLV`).

---

## Phase 3: User Story 1 - Full Physical Radiative Solver Integration (Priority: P1) 🎯 MVP

**Goal**: Port OPMIE and replace the mock flux loop inside the engine orchestrator.

- [x] T003 Port the `OPMIE` light propagation and layer integration subroutine from `cldj_fjx_sub_mod.F90` to `include/cloudj/radiative_solver.hpp`.
- [x] T004 Wire up the ported `OPMIE` and `MIESCT` routines directly inside `CloudJ::Engine::calculate_photolysis_rates` in `include/cloudj/cloudj.hpp`, completely replacing the mock direct flux scaling loop.

---

## Phase 4: User Story 2 - High-Order Polynomial Exponential Optimization (Priority: P1)

**Goal**: Implement high-order Minimax polynomial exponential approximations.

- [x] T005 Implement a 4th-degree Minimax polynomial approximation to evaluate `std::exp(-tau)` beam extinctions inside `OPMIE` in `include/cloudj/radiative_solver.hpp`.

---

## Phase 5: Polish & Cross-Cutting Concerns

**Purpose**: Validation and benchmarking.

- [x] T006 [TDD] Compile and run unit tests (`test_library_api.cpp`) to verify API changes and numerical parity.
- [x] T007 Run the 500,000-run property fuzzer to guarantee $10^{-12}$ exact numerical parity remains intact.
- [x] T008 Run the high-resolution benchmarking script (`run_benchmarks.py`) to verify the speedups.

---

## Dependencies & Execution Order

1. **Phase 1 (Setup)**: N/A
2. **Phase 2 (Foundational)**: Port `LEGND0` and `MIESCT`.
3. **Phase 3 (User Story 1)**: Port `OPMIE` and wire up inside `Engine`.
4. **Phase 4 (User Story 2)**: Minimax exponential optimization inside `OPMIE`.
5. **Phase 5 (Polish)**: Final extensive validation.
