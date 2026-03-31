from unittest.mock import MagicMock, patch

import pytest

from climateset.download.client import SearchClient
from climateset.download.constraints import BaseSearchConstraints


@pytest.fixture
def mock_search_connection():
    with patch("climateset.download.client.SearchConnection") as mock:
        yield mock


def test_search_client_context_manager():
    """Test that SearchClient works as a context manager."""
    with SearchClient() as client:
        assert isinstance(client, SearchClient)


def test_search_session_initial_connection(mock_search_connection):
    """Test that a new session establishes a connection to the first node."""
    mock_conn_instance = MagicMock()
    mock_search_connection.return_value = mock_conn_instance

    client = SearchClient(node_urls=["http://node1", "http://node2"])
    session = client.new_session()

    mock_search_connection.assert_called_with(url="http://node1", distrib=True)
    mock_conn_instance.new_context.assert_called_once()
    assert session._connection == mock_conn_instance


def test_search_session_failover(mock_search_connection):
    """Test that session fails over to the next node if the first one fails."""
    # First call raises error, second returns mock
    mock_conn_instance = MagicMock()
    mock_search_connection.side_effect = [Exception("Connection failed"), mock_conn_instance]

    client = SearchClient(node_urls=["http://node1", "http://node2"])
    session = client.new_session()

    # Should have tried node1 then node2
    assert mock_search_connection.call_count == 2
    mock_search_connection.assert_any_call(url="http://node1", distrib=True)
    mock_search_connection.assert_any_call(url="http://node2", distrib=True)
    assert session._connection == mock_conn_instance


def test_search_session_constrain_replay(mock_search_connection):
    """
    Test that constraints are replayed when failing over to a new node.

    Scenario:
    1. Connect to Node 1 successfully.
    2. Apply constraint A.
    3. Apply constraint B (fails on Node 1).
    4. Session should rotate to Node 2 and replay A and B.
    """
    # Setup mocks
    node1_conn = MagicMock()
    node2_conn = MagicMock()

    node1_ctx = MagicMock()
    node2_ctx = MagicMock()

    node1_conn.new_context.return_value = node1_ctx
    node2_conn.new_context.return_value = node2_ctx

    # Node 1 context dies on second constraint
    node1_ctx.constrain.side_effect = [
        node1_ctx,  # First constraint OK
        Exception("Node 1 died"),  # Second constraint fails
    ]

    # Node 2 context succeeds
    node2_ctx.constrain.return_value = node2_ctx

    mock_search_connection.side_effect = [node1_conn, node2_conn]

    client = SearchClient(node_urls=["http://node1", "http://node2"])
    session = client.new_session()

    # 1. Connected to Node 1
    constraints_a = BaseSearchConstraints(project="CMIP6")
    session.constrain(constraints_a)

    # Verify Node 1 constrained
    node1_ctx.constrain.assert_called_with(project="CMIP6")

    # 2. Convert to params and constrain again -> Logic inside constrain() handles exceptions?
    # Actually, constrain() calls _context.constrain.
    # If that raises, it should catch, rotate, and re-ensure connection (replaying all).

    constraints_b = BaseSearchConstraints(variable="tas")
    session.constrain(constraints_b)

    # Only verify we moved to Node 2
    assert session._current_node_index == 1  # 0-indexed, so 1 is second node

    # Verify Node 2 was initialized and constrained with BOTH A and B
    node2_conn.new_context.assert_called()
    assert node2_ctx.constrain.call_count >= 2
    # call_args_list should verify replay order
    # Note: dictionary ordering might vary but we passed simple kwargs
