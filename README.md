# AERO740 Quadcopter Project Collaboration

## Project Overview

This project investigates the **control and behavior of a quadcopter operating in the vertical plane** under realistic failure modes and complications. The research will introduce complexity of the following types:

1. **Communication Time Delay** – Latency between controller computation and actuator execution
2. **Unobservable Disturbance** – Vertical wind profile modeled as an external disturbance
3. **Model Mismatch** – Difference between controller prediction dynamics and true vehicle performance
4. **State and Input Constraints** – Obstacle avoidance, position limits, and time-of-arrival requirements

### Objectives

The goal is to **explore how these failure modes manifest within a Model Predictive Control (MPC) framework** and to evaluate mitigation strategies that preserve both feasibility and performance.

### Mathematical Framework

The quadcopter dynamics are approximated as a **linear continuous-time system**:

$$\dot{x}(t) = A x(t) + B u(t) + w_c(t)$$

where $w_c(t)$ represents external wind disturbance.

The MPC controller solves a finite-horizon optimal control problem over prediction horizon $\Delta T_H$, discretized into $N$ steps ($\Delta T_k = \Delta T_H / N$), using the **discrete-time prediction model**:

$$x_{k+1} = A_d x_k + B_d u_k + w_{d,k}$$

This framework enables systematic investigation of how different complications affect control performance and feasibility.

---

## Getting Started with Git

Our team uses Git for version control. If you're new to Git, here's how to get started with this project.

### 1. Clone the Repository

**For the first time, clone the repository to your local machine:**

```bash
git clone https://github.com/omokeefe-umich/aero740_project_colab.git
cd aero740_project_colab
```

**Learn more:** [Git Clone Documentation](https://git-scm.com/docs/git-clone)

### 2. Check the Status

**Before making changes, always check the current status:**

```bash
git status
```

This shows which files have been modified, added, or deleted.

**Learn more:** [Git Status Documentation](https://git-scm.com/docs/git-status)

### 3. Pull Latest Changes

**Always pull the latest changes before starting work:**

```bash
git pull origin main
```

This ensures you have the most recent work from your teammates and prevents merge conflicts.

**Learn more:** [Git Pull Documentation](https://git-scm.com/docs/git-pull)

### 4. Create a Branch (Recommended)

**For larger feature work, create a feature branch:**

```bash
git checkout -b feature/your-feature-name
```

This keeps work organized and allows teammates to work independently.

**Learn more:** [Git Checkout Documentation](https://git-scm.com/docs/git-checkout)

### 5. Make Changes and Stage Them

**After modifying files, stage your changes:**

```bash
git add .
```

Or stage specific files:

```bash
git add path/to/file.m
```

**Learn more:** [Git Add Documentation](https://git-scm.com/docs/git-add)

### 6. Commit Your Changes

**Write a clear, descriptive commit message:**

```bash
git commit -m "Add description of what you changed"
```

Good commit messages help teammates understand the history of changes.

**Learn more:** [Git Commit Documentation](https://git-scm.com/docs/git-commit)

### 7. Push Your Changes

**Push your commits to the remote repository:**

```bash
git push origin feature/your-feature-name
```

Or to the main branch:

```bash
git push origin main
```

**Learn more:** [Git Push Documentation](https://git-scm.com/docs/git-push)

### 8. Create a Pull Request

**When your feature is complete, create a pull request for review** (if using GitHub/GitLab/Bitbucket).

This allows teammates to review your changes before merging.

---

## Common Git Workflows

### Daily Workflow

```bash
# Start your day
git pull origin main

# Make changes to your files
# ...

# Check what you changed
git status

# Stage and commit
git add .
git commit -m "Clear description of changes"

# Push to remote
git push origin main
```

### Undoing Changes

**Discard local changes to a file:**

```bash
git checkout -- path/to/file.m
```

**Undo the last commit (but keep changes):**

```bash
git reset --soft HEAD~1
```

**View commit history:**

```bash
git log --oneline
```

**Learn more:** [Git Reset Documentation](https://git-scm.com/docs/git-reset)

---

## Helpful Resources

- **[Git Official Documentation](https://git-scm.com/doc)** – Comprehensive Git reference
- **[Git Cheat Sheet](https://github.github.com/training-kit/downloads/github-git-cheat-sheet.pdf)** – Quick reference guide
- **[GitHub Learning Lab](https://lab.github.com/)** – Interactive Git tutorials
- **[Atlassian Git Tutorials](https://www.atlassian.com/git/tutorials)** – Practical Git guides
- **[Pro Git Book](https://git-scm.com/book/en/v2)** – Free comprehensive Git guide

---

## Tips for Collaboration

1. **Pull before you push** – Always fetch the latest changes before pushing
2. **Write clear commit messages** – Help your teammates understand your changes
3. **Commit often** – Small, logical commits are easier to understand and revert if needed
4. **Communicate with your team** – Let teammates know what you're working on
5. **Review changes before committing** – Use `git diff` to see exactly what you changed

---

## Questions?

If you're stuck or have questions about Git:
- Check the [Git Official Documentation](https://git-scm.com/doc)
- Ask your teammates
- Search Stack Overflow for common issues
