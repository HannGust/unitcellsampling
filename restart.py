from persistqueue import Queue


class RestartQueue:
    def __init__(self, name='restartgrid'):
        self.name = name
        self.queue = Queue(name, autosave=True)

    def is_empty(self):
        return self.queue.info['size'] == 0

    def get(self):
        return self.queue.get()

    def put(self, data):
        return self.queue.put(data)
