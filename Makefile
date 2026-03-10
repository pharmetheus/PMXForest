.PHONY: check-lint
check-lint: export LINTR_ERROR_ON_LINT = true
check-lint:
	R --slave --no-save --no-restore -e 'lintr::lint_package()'

.PHONY: check-format
check-format:
	R --slave --no-save --no-restore -e 'styler::style_pkg()'
	@git --no-pager diff --color=always -- . > styler.diff || true
	if [ -s styler.diff ]; then \
	  echo "The following files need styling:"; \
	  cat styler.diff; \
	  git checkout -- .; \
	  exit 1; \
	else \
	  echo "Styling OK"; \
	fi

.PHONY: test
test:
	R --slave --no-save --no-restore -e 'testthat::set_max_fails(Inf); testthat::test_local(stop_on_failure=TRUE)'
