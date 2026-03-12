# ACP Stdio Mode

DrugClaw can run as an Agent Client Protocol (ACP) server over stdio:

```sh
drugclaw acp
```

## When to use it

Use ACP mode when another local tool wants to treat DrugClaw as a sessioned chat runtime over stdio instead of using a chat adapter or the Web API.

Typical cases:

- local editor or IDE integrations
- terminal wrappers that want ACP transport
- local automation that already speaks ACP

## Behavior

- uses the normal `drugclaw.config.yaml`
- persists ACP conversations through the standard runtime storage
- supports `/stop` to cancel the active run for the ACP session
- keeps the normal tool loop and provider stack

## Verification

1. Run `drugclaw doctor`.
2. Start `drugclaw acp`.
3. Connect with an ACP client.
4. Send one prompt and confirm a normal response.
5. Send a follow-up prompt in the same session and confirm context is preserved.
6. Trigger a long-running request, then send `/stop` and confirm cancellation works.

## Related docs

- `README.md`
- `TEST.md`
