#!/usr/bin/env bash

set -Eeuo pipefail

readonly app_url="http://127.0.0.1:3838/"
readonly server_log="/tmp/targetscore-shiny.log"
readonly response_file="/tmp/targetscore-shiny-response.html"
readonly pipeline_output="/tmp/targetscore-vignette-output"
server_pid=""

print_server_log() {
    if [[ -f "${server_log}" ]]; then
        printf '%s\n' "--- Shiny server log ---" >&2
        cat "${server_log}" >&2
    fi
}

stop_server() {
    if [[ -z "${server_pid}" ]]; then
        return
    fi

    if ! kill -0 "${server_pid}" 2>/dev/null; then
        wait "${server_pid}" 2>/dev/null || true
        server_pid=""
        return
    fi

    kill -TERM "${server_pid}" 2>/dev/null || true
    for _ in $(seq 1 50); do
        if ! kill -0 "${server_pid}" 2>/dev/null; then
            wait "${server_pid}" || true
            server_pid=""
            return
        fi
        sleep 0.1
    done

    kill -KILL "${server_pid}" 2>/dev/null || true
    wait "${server_pid}" 2>/dev/null || true
    server_pid=""
}

cleanup() {
    local exit_status=$?
    if ((exit_status != 0)); then
        print_server_log
    fi
    stop_server
    exit "${exit_status}"
}

trap cleanup EXIT

Rscript -e "shiny::runApp(system.file('shiny', package = 'targetscore'), host = '0.0.0.0', port = 3838)" \
    >"${server_log}" 2>&1 &
server_pid=$!

server_ready=false
for _ in $(seq 1 60); do
    if ! kill -0 "${server_pid}" 2>/dev/null; then
        wait "${server_pid}" || server_status=$?
        server_pid=""
        printf 'Shiny server exited before readiness (status %s).\n' "${server_status:-0}" >&2
        exit 1
    fi

    if curl --fail --silent --max-time 2 "${app_url}" >"${response_file}"; then
        server_ready=true
        break
    fi
    sleep 1
done

if [[ "${server_ready}" != true ]]; then
    printf 'Shiny server was not ready at %s within 60 seconds.\n' "${app_url}" >&2
    exit 1
fi

if ! grep --quiet "Target Score" "${response_file}"; then
    printf 'Shiny root endpoint did not contain the expected application title.\n' >&2
    exit 1
fi

stop_server

mkdir -p "${pipeline_output}"
Rscript -e "rmarkdown::render('/opt/targetscore/vignettes/target_score_tutorial.Rmd', output_format = 'html_document', output_dir = '${pipeline_output}', params = list(output_dir = '${pipeline_output}'))"

if [[ ! -s "${pipeline_output}/target_score_tutorial.html" ]]; then
    printf 'Vignette pipeline did not produce target_score_tutorial.html.\n' >&2
    exit 1
fi
