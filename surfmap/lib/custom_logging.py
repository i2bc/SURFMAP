class  Log:
    def info(self, message):
        print(message)

    def section(self, message):
        print(f"\n\n{message}")
        print(f"{len(message)*'-'}\n")
