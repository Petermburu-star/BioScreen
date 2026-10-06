"""Launch BioScreen Studio from a Jupyter notebook."""
import threading
import time

from IPython.display import IFrame

from bioscreen.studio.app import app, _load_engine

_server_started = False


def launch(port=5000, height=1100):
    """Start the Flask server in a background thread and embed it in the notebook."""
    global _server_started

    def _run():
        _load_engine()
        app.run(host="127.0.0.1", port=port, debug=False, use_reloader=False)

    if not _server_started:
        threading.Thread(target=_run, daemon=True).start()
        _server_started = True
        time.sleep(2)

    return IFrame(f"http://127.0.0.1:{port}", width="100%", height=height)
