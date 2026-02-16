from typing import Any, List, Optional

from pyesgf.search import SearchConnection
from pyesgf.search.context import DatasetSearchContext

from climateset.download.constants import NODE_LINK_URLS
from climateset.download.constraints import BaseSearchConstraints
from climateset.utils import create_logger

LOGGER = create_logger(__name__)


class SearchClient:
    """
    Client for performing searches against ESGF nodes with failover support.

    Acts as a factory for SearchSession objects.
    """

    def __init__(self, node_urls: List[str] | None = None, distrib: bool = True):
        self.node_urls = node_urls if node_urls is not None else NODE_LINK_URLS
        self.distrib = distrib
        self.logger = LOGGER

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def new_session(self) -> "SearchSession":
        """Start a new search session."""
        return SearchSession(self.node_urls, self.distrib, self.logger)


class SearchSession:
    """
    Stateful session for building iterative search queries.

    Handles node failover by replaying applied constraints.
    """

    def __init__(self, node_urls: List[str], distrib: bool, logger):
        self.node_urls = node_urls
        self.distrib = distrib
        self.logger = logger

        # History of constraints applied to this session
        self._constraints_history: List[BaseSearchConstraints] = []

        # State relative to correct active connection
        self._current_node_index = 0
        self._connection: Optional[SearchConnection] = None
        self._context: Optional[DatasetSearchContext] = None

        # Initialize connection logic
        self._ensure_connection()

    def _ensure_connection(self):
        """
        Ensures a valid connection/context exists.

        If not, attempts to connect to available nodes. Once connected, replays history.
        """
        if self._context is not None:
            return

        while self._current_node_index < len(self.node_urls):
            url = self.node_urls[self._current_node_index]
            try:
                self.logger.info(f"Connecting to ESGF node: {url}")
                self._connection = SearchConnection(url=url, distrib=self.distrib)

                # Create fresh context
                ctx = self._connection.new_context()

                # Replay constraints
                for constraints in self._constraints_history:
                    params = constraints.to_esgf_params()
                    if params:
                        ctx = ctx.constrain(**params)

                self._context = ctx
                return
            except Exception as e:  # pylint: disable=broad-exception-caught
                self.logger.warning(f"Failed to connect to {url}: {e}")
                self._current_node_index += 1
                self._connection = None
                self._context = None

        raise ConnectionError(f"Could not connect to any ESGF node. Tried: {self.node_urls}")

    def _rotate_node(self):
        """Force rotation to the next node (e.g. after a search failure)."""
        self.logger.info("Rotating to next ESGF node...")
        self._current_node_index += 1
        self._connection = None
        self._context = None
        self._ensure_connection()

    def constrain(self, constraints: BaseSearchConstraints) -> "SearchSession":
        """Apply a new set of constraints to the session."""
        self._constraints_history.append(constraints)

        # If we have an active context, apply immediately.
        # If not (e.g. all nodes down), _ensure_connection will handle it next time.
        if self._context:
            params = constraints.to_esgf_params()
            if params:
                try:
                    self._context = self._context.constrain(**params)
                except Exception as e:  # pylint: disable=broad-exception-caught
                    self.logger.warning(f"Error applying constraints on current node: {e}")
                    self._rotate_node()
        else:
            # Try to establish connection if we were disconnected
            try:
                self._ensure_connection()
            except ConnectionError:
                pass  # Delay error until actual search/facet request

        return self

    def get_available_facets(self, facet_name: str) -> List[str]:
        """
        Get available counts/values for a specific facet.

        Retries on other nodes if current fails.
        """
        max_attempts = len(self.node_urls)
        attempts = 0

        while attempts < max_attempts:
            try:
                self._ensure_connection()
                if facet_name in self._context.facet_counts:
                    return list(self._context.facet_counts[facet_name].keys())
                return []
            except Exception as e:  # pylint: disable=broad-exception-caught
                self.logger.warning(f"Error fetching facets from {self.node_urls[self._current_node_index]}: {e}")
                self._rotate_node()
                attempts += 1

        return []

    def search(self) -> List[Any]:
        """
        Execute the search using applied constraints.

        Retries on other nodes if current fails.
        """
        max_attempts = len(self.node_urls)
        attempts = 0

        while attempts < max_attempts:
            try:
                self._ensure_connection()
                return self._context.search()
            except Exception as e:  # pylint: disable=broad-exception-caught
                self.logger.warning(f"Search failed on {self.node_urls[self._current_node_index]}: {e}")
                self._rotate_node()
                attempts += 1

        raise ConnectionError("Search failed on all available nodes.")
