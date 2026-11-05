CXX ?= g++
CC ?= gcc
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pthread -fopenmp
CPPFLAGS += -Isrc -Ithird_party/minimap2 -Ithird_party/htsbox/htslib
LDLIBS += -lz -lm -pthread
BUILD := build

MM2 := third_party/minimap2/libminimap2.a
HTS := third_party/htsbox/htslib/libhts.a
SOURCES := main Inference InputGraph PairedReadBinning DipcallRunner PanGenieVcf gfa-base AssemblyCalling VariantOutput
OBJECTS := $(addprefix $(BUILD)/,$(addsuffix .o,$(SOURCES)))
CALLER_OBJECTS := $(BUILD)/EmbeddedPileup.o $(addprefix $(BUILD)/htsbox-,faidx.o razf.o bedidx.o)
CALLER_LIBS := $(BUILD)/AssemblyCalling.o $(BUILD)/VariantOutput.o $(CALLER_OBJECTS) $(MM2) $(HTS)

TESTS := test wide_test test_sequence_index \
	test_parallel_dp test_paired_bins test_constants test_global_phase \
	test_dense_pair_repair test_locus_dosage test_library_roles \
	test_target_counts test_observed_dosage test_complete_dosage test_window_selection test_copy_levels \
	test_free_swaps test_het_vote test_native_localizer test_native_allelic \
	test_native_aliases test_native_graph_localizer test_native_read_evidence test_native_hic_split \
	test_native_factors test_graph_loader test_sharded_words test_hot_containers test_factor_batches \
	test_exact_match test_progress_log
TEST_BINS := $(addprefix $(BUILD)/,$(TESTS))
TOY := python3 test/joint/raw_toy_test.py --executable ./DipGenie --recombinant

.PHONY: all check clean
all: DipGenie

$(BUILD):
	mkdir -p $(BUILD)
$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
$(BUILD)/EmbeddedPileup.o: src/EmbeddedPileup.c third_party/htsbox/pileup.c | $(BUILD)
	$(CC) -O3 -std=gnu11 -Wall -Wextra -Ithird_party/htsbox/htslib -MMD -MP -c $< -o $@
$(BUILD)/htsbox-%.o: third_party/htsbox/%.c | $(BUILD)
	$(CC) -O3 -std=gnu11 -Wall -Ithird_party/htsbox/htslib -MMD -MP -c $< -o $@
$(MM2):
	$(MAKE) -C third_party/minimap2 libminimap2.a CFLAGS='-Wall -O3 -Wc++-compat -Dks_init=dg_mm_ks_init -Dks_destroy=dg_mm_ks_destroy -Dks_getuntil2=dg_mm_ks_getuntil2'
$(HTS):
	$(MAKE) -C third_party/htsbox/htslib libhts.a

DipGenie: $(OBJECTS) $(CALLER_OBJECTS) $(MM2) $(HTS)
	$(CXX) $(CXXFLAGS) $^ $(LDLIBS) -o $@

# Tests
$(TEST_BINS): $(wildcard third_party/unordered_dense/ankerl/*.h)
$(BUILD)/test_assembly_calling: test/joint/test_assembly_calling.cpp $(CALLER_LIBS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(CALLER_LIBS) $(LDLIBS) -o $@
$(BUILD)/test_paired_bins: test/joint/test_paired_bins.cpp $(BUILD)/PairedReadBinning.o | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/PairedReadBinning.o $(LDLIBS) -o $@
$(BUILD)/call_fasta: data/expts/Chromosome_1/05_baselines/call_fasta.cpp $(BUILD)/DipcallRunner.o $(CALLER_LIBS) $(wildcard src/*.hpp) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/DipcallRunner.o $(CALLER_LIBS) $(LDLIBS) -o $@
$(BUILD)/%: test/joint/%.cpp $(wildcard src/*.hpp) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDLIBS) -o $@

check: DipGenie $(TEST_BINS) $(BUILD)/test_assembly_calling $(BUILD)/call_fasta
	@set -e; for t in $(TEST_BINS); do "$$t"; done
	@set -e; for lib in sr hic mixed; do for graph in "" --cyclic; do \
		$(TOY) --libraries $$lib $$graph; \
		$(TOY) --libraries $$lib --unobserved-alleles $$graph; \
		if test $$lib != sr; then $(TOY) --libraries $$lib --cross-locus $$graph; fi; \
	done; done
	$(TOY) --libraries mixed --dipcall
	python3 test/joint/test_inputs.py ./DipGenie
	python3 test/joint/assembly_parity.py $(BUILD)/test_assembly_calling --fasta-helper $(BUILD)/call_fasta
	python3 test/joint/variant_output_parity.py $(BUILD)/test_assembly_calling
	python3 data/expts/Chromosome_1/04_dipgenie/test_run.py
	python3 data/expts/Chromosome_1/06_evaluation/test_evaluate.py
	python3 data/expts/Chromosome_1/06_evaluation/test_locked_cache.py

clean:
	rm -rf $(BUILD) DipGenie

-include $(OBJECTS:.o=.d) $(CALLER_OBJECTS:.o=.d)
